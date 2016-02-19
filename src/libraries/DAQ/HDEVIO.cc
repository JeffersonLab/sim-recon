// $Id$
//
//    File: HDEVIO.cc
// Created: Wed Dec 10 07:22:00 EST 2014
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#include <stdlib.h>
#include <string.h>

#include "HDEVIO.h"

//---------------------------------
// HDEVIO    (Constructor)
//---------------------------------
HDEVIO::HDEVIO(string filename):filename(filename)
{
	is_open = false;
	ifs.open(filename.c_str());
	if(!ifs.is_open()){
		ClearErrorMessage();
		err_mess << "Unable to open EVIO file: " << filename;
		return;
	}
	
	fbuff_size = 10000000; // 40MB input buffer
	fbuff = new uint32_t[fbuff_size];;
	fnext = fbuff;
	fbuff_end = fbuff;
	fbuff_len = 0;
	_gcount = 0;
	
	buff_limit = 5000000; // Don't allow us to allocate more than 5M words for read buffer
	buff_size = 1024;     // initialize with 4kB buffer
	buff_len = 0;
	buff = new uint32_t[buff_size];
	next = buff; // needed so initial calculation of left is 0
	last_event_len = 0;
	err_code = HDEVIO_OK;

	Nblocks = 0;
	Nevents = 0;
	Nerrors = 0;
	Nbad_blocks = 0;
	Nbad_events = 0;
	
	is_open = true;
}

//---------------------------------
// buff_read
//---------------------------------
void HDEVIO::buff_read(char* s, streamsize nwords)
{
	/// Read n bytes into the specified buffer.
	/// This will read from the file if needed, or simply copy
	/// the bytes from the existing fbuff. It is purposely meant
	/// as a drop-in replacement for ifstream::read so that we
	/// can make actual reads from the file in larger chunks and
	/// avoid backwards seeks on the actual file. This is needed
	/// because the Lustre filesystem is optimized for large data
	/// reads but smaller reads with backwards seeks tend to cripple
	/// its performance. The main difference between this and
	/// ifstream::read is that this does not return an istream&
	/// reference.
	
	// Number of words left in fbuff	
	uint64_t left = ((uint64_t)fbuff_end - (uint64_t)fnext)/sizeof(uint32_t);
	
	// First, copy what's already in fbuff
	uint64_t Ncopied = nwords<(int64_t)left ? (uint64_t)nwords:left;
	_gcount = Ncopied*sizeof(uint32_t);
	if(_gcount>0) memcpy((char*)s, (char*)fnext, _gcount);
	left -= Ncopied;
	fnext += Ncopied;
	s += _gcount; // advance pointer to user buff in case we need to write more
	
	// If needed, read in another chunk from the file.
	// Try and keep last 8 words so if a seekg is called to back
	// us up that amount, we don't have to call seekg on the file.
	if( (int64_t)Ncopied < nwords ){
	
		// Initialize to start of buffer
		uint32_t *myfbuff = fbuff;
		uint32_t myfbuff_size = fbuff_size;
		fbuff_len = 0;
	
		// If at least 8 words exist in fbuff, then copy last
		// 8 words into front of fbuff.
		if( fbuff_len >= 8 ){
			memcpy((char*)fbuff, (char*)&fbuff[fbuff_len-8], 8*sizeof(uint32_t));
			myfbuff = &fbuff[8];
			myfbuff_size -= 8;
			fbuff_len += 8;
		}

		// Read in chunk from file
		ifs.read((char*)myfbuff, myfbuff_size*sizeof(uint32_t));
		fbuff_len += (uint64_t)(ifs.gcount()/sizeof(uint32_t));
		fnext = myfbuff;
		fbuff_end = &fbuff[fbuff_len];

		// Copy remainder of request
		uint64_t myleft = ((uint64_t)fbuff_end - (uint64_t)fnext)/sizeof(uint32_t);
		
		uint64_t mynwords = nwords - Ncopied;
		uint64_t myNcopied = mynwords<myleft ? mynwords:myleft;
		uint64_t mygcount = myNcopied*sizeof(uint32_t);
		if(mygcount>0) memcpy((char*)s, (char*)fnext, mygcount);
		fnext += myNcopied;
		_gcount += myNcopied*sizeof(uint32_t);		
	}
}

//---------------------------------
// buff_seekg
//---------------------------------
void HDEVIO::buff_seekg (streamoff off, ios_base::seekdir way)
{
	// Convert offset from bytes to words
	int64_t off_words = (int64_t)off/(int64_t)sizeof(uint32_t);

	// find current position relative to start of buffer in words
	int64_t fpos = (int64_t)((uint64_t)fnext - (uint64_t)fbuff)/sizeof(uint32_t);

	// Seek depending on what user has requested
	if( way == ios_base::cur ){
	
		// Add requested offset
		fpos += off_words; // desired position relative to start of fbuff
		
		if(fpos>=0 && fpos<(int64_t)fbuff_len){
			// Request is to a point inside current buffer
			fnext = &fbuff[fpos];
			_gcount = 0;
		}else if(fpos<0){
			// Seek point is outside of buffer. Move file pointer
			// and indicate buffer is now empty so next buff_read()
			// call will force a read.
			
			// Current file position should be just after end of fbuff
			// Subtract fbuff_len to move it back to start of fbuff.
			// fpos is position relative to start of fbuff.
			off = (streamoff)fpos - (streamoff)fbuff_len; // offset relative to actual current file position
			
			ifs.seekg(off);
			
			// Set fbuff parameters to indicate no valid data
			fnext = fbuff;
			fbuff_end = fbuff;
			fbuff_len = 0;
			_gcount = 0;

		}
		
	}else{
		_DBG_<<"buff_seekg called with something other than ios_base::cur is unsupported!" << endl;
		exit(-1);
	}
}

//---------------------------------
// ReadBlock
//---------------------------------
bool HDEVIO::ReadBlock(void)
{
	/// Read in the next EVIO block. Return true if successful
	/// and false otherwise.
		
	err_code = HDEVIO_OK;
	next = NULL;

	buff_read((char*)buff, 8);
	uint32_t valid_words = buff_gcount()/sizeof(uint32_t);
	if(valid_words != 8){
		SetErrorMessage("Could not read in 8 word EVIO block header!");
		err_code = HDEVIO_FILE_TRUNCATED;
		Nerrors++;
		return false;
	}

	// Check endianess
	if(buff[7]!=0xc0da0100 && buff[7]!=0x0001dac0){
		ClearErrorMessage();
		err_mess << "Magic word not valid!: " << HexStr(buff[7]) << endl;
		err_code = HDEVIO_BAD_BLOCK_HEADER;
		Nerrors++;
		Nbad_blocks++;
		return false;
	}
	swap_needed = (buff[7]==0x0001dac0);

	// Swap header (if needed)
	if(swap_needed) swap_block(buff, 8, buff);

	// Re-allocate buffer if needed so we can read in entire block
	uint32_t block_length = buff[0];
	if(buff_size < block_length){
		if(block_length > buff_limit){
			ClearErrorMessage();
			err_mess << "ERROR: EVIO block length greater than allocation limit (" << block_length <<" > " << block_length << " words)" << endl;
			err_code = HDEVIO_BLOCKSIZE_GREATER_THAN_LIMIT;
			Nerrors++;
			return false;
		}
		if(buff) delete[] buff;
		buff_size = block_length;
		buff = new uint32_t[buff_size];
		if(buff == NULL){
			ClearErrorMessage();
			err_mess << "ERROR: unable to allocate " << block_length <<" words" << endl;
			err_code = HDEVIO_MEMORY_ALLOCATION_ERROR;
			return false;
		}
		
		// Re-read in the block header
		buff_seekg(-8*sizeof(uint32_t), ifs.cur);
		buff_read((char*)buff, 8);
		if(swap_needed) swap_block(buff, 8, buff);
	}
	
	if(block_length == 8){
		// block_length =8 indicates end of file.
		SetErrorMessage("end of file");
		err_code = HDEVIO_EOF;
		return false;
	}

	// Read payload of block
	buff_read((char*)&buff[8], (block_length-8));
	valid_words = 8 + buff_gcount()/sizeof(uint32_t);
	if(valid_words < block_length){
		ClearErrorMessage();
		err_mess << "Error reading in EVIO entire block! (block number: " << buff[1] << ")" << endl;
		err_mess << "valid_words="<<valid_words << " block_length=" << block_length;
		err_code = HDEVIO_FILE_TRUNCATED;
		Nerrors++;
		return false;
	}
	
	// Set pointers
	buff_len = valid_words;
	buff_end = &buff[valid_words];
	next = &buff[8];
	bh = (BLOCKHEADER_t*)buff;
	
	Nblocks++;
	return true;
}

//---------------------------------
// read
//---------------------------------
bool HDEVIO::read(uint32_t *user_buff, uint32_t user_buff_len)
{
	/// Read the next EVIO event into the user supplied buffer.
	/// Return true if successful and false otherwise. Details of
	/// the error will be in err_mess.
	
	err_code = HDEVIO_OK;

	// calculate remaining valid words in buffer
	uint32_t left = buff_len - (uint32_t)(((unsigned long)next - (unsigned long)buff)/sizeof(uint32_t));

	// Read in another event block if necessary
	if(left < 1 || next==NULL){
		bool isgood = ReadBlock();
		if(!isgood) return false;
		left = buff_len - 8;
	}
	
	if(next == NULL){
		SetErrorMessage("No valid events in buffer");
		err_code = HDEVIO_NO_EVENTS_IN_BUFFER;
		return false;
	}
	
	// Check if next event will fit into user supplied buffer
	uint32_t event_len = next[0];
	if(swap_needed) swap_block(&event_len, 1, &event_len);
	event_len++; // include length word for EVIO bank
	last_event_len = event_len;

	// Check that event isn't claiming to be larger than EVIO block
	if( event_len > left ){
		ClearErrorMessage();
		err_mess << "WARNING: EVIO bank indicates a bigger size than block header (" << event_len << " > " << left << ")";
		next = &buff[buff_len]; // setup so subsequent call will read in another block
		err_code = HDEVIO_EVENT_BIGGER_THAN_BLOCK;
		Nerrors++;
		Nbad_blocks++;
		return false;
	}

	// Check if user buffer is big enough to hold this 
	if(event_len > user_buff_len){
		ClearErrorMessage();
		err_mess << "user buffer too small for event (" << user_buff_len << " < " << event_len << ")";
		err_code = HDEVIO_USER_BUFFER_TOO_SMALL;
		return false;
	}
	
	// Copy entire event into user buffer, swapping if needed during copy
	bool isgood = true;
	if(swap_needed){
		uint32_t Nswapped = swap_bank(user_buff, next, event_len);
		isgood = (Nswapped == event_len);
	}else{
		memcpy(user_buff, next, event_len*sizeof(uint32_t));
	}
	
	// Advance next pointer to next EVIO event or past end of buffer.
	next = &next[event_len];
	
	if(isgood) Nevents++;

	return isgood;
}


//---------------------------------
// ~HDEVIO    (Destructor)
//---------------------------------
HDEVIO::~HDEVIO()
{
	if(ifs.is_open()) ifs.close();
	if(buff) delete[] buff;
}

//---------------------------------
// swap_bank
//---------------------------------
uint32_t HDEVIO::swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO bank. If the bank contains data, it is automatically
	/// swapped according to it's type. If the bank is a container of other
	/// containers, then this repeatedly calls the swapper methods for the
	/// appropriate container type (bank, tagsegment, segment). This means
	/// that this method will be recursive in the cases where it is a bank
	/// of banks.

	if(len < 2){
		SetErrorMessage("Attempt to swap bank with len<2");
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	// Swap length and header words
	swap_block(inbuff, 2, outbuff);
	uint32_t bank_len = outbuff[0];
	if((bank_len+1) > len){
		ClearErrorMessage();
		err_mess << "WARNING: Bank length word exceeds valid words in buffer (" << bank_len+1 << " > " << len << ")";
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	uint32_t type = (outbuff[1]>>8) & 0xFF;
	uint32_t Nwords = bank_len - 1;
	uint32_t Nswapped = 2;
	switch(type){
		case 0x0a:  // 64 bit unsigned int
		case 0x08:  // 64 bit double
		case 0x09:  // 64 bit signed int
			swap_block((uint64_t*)&inbuff[2], Nwords/2, (uint64_t*)&outbuff[2]);
			Nswapped += Nwords;
			break;
		case 0x01:  // 32 bit unsigned int
		case 0x02:  // 32 bit float
		case 0x0b:  // 32 bit signed int
			swap_block(&inbuff[2], Nwords, &outbuff[2]);
			Nswapped += Nwords;
			break;
		case 0x05:  // 16 bit unsigned int
		case 0x04:  // 16 bit signed int
			swap_block((uint16_t*)&inbuff[2], Nwords*2, (uint16_t*)&outbuff[2]);
			Nswapped += Nwords;
			break;
		case 0x00:  // 32 bit unknown (not swapped)
		case 0x07:  // 8 bit unsigned int
		case 0x06:  // 8 bit signed int
			memcpy((uint8_t*)&outbuff[2], (uint8_t*)&inbuff[2], Nwords*sizeof(uint32_t));
			Nswapped += Nwords;
			break;
		case 0x0c:
			while(Nswapped < (Nwords+2)){
				uint32_t N = swap_tagsegment(&outbuff[Nswapped], &inbuff[Nswapped], (Nwords+2)-Nswapped);
				if(N == 0) return Nswapped;
				Nswapped += N;
			}
			break;
		case 0x0d:
		case 0x20:
			while(Nswapped < (Nwords+2)){
				uint32_t N = swap_segment(&outbuff[Nswapped], &inbuff[Nswapped], (Nwords+2)-Nswapped);
				if(N == 0) return Nswapped;
				Nswapped += N;
			}
			break;
		case 0x0e:
		case 0x10:
			while(Nswapped < (Nwords+2)){
				uint32_t N = swap_bank(&outbuff[Nswapped], &inbuff[Nswapped], (Nwords+2)-Nswapped);
				if(N == 0) return Nswapped;
				Nswapped += N;
			}
			break;
		default:
			ClearErrorMessage();
			err_mess << "WARNING: unknown bank type (0x" << hex << type << dec << ")";
			err_code = HDEVIO_UNKNOWN_BANK_TYPE;
			Nerrors++;
			Nbad_events++;
			return 0;
			break;
	}

	return Nswapped;
}

//---------------------------------
// swap_tagsegment
//---------------------------------
uint32_t HDEVIO::swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO tagsegment. 

	if(len < 1){
		SetErrorMessage("Attempt to swap segment with len<1");
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	// Swap header/length word
	swap_block(inbuff, 1, outbuff);
	uint32_t bank_len = outbuff[0] & 0xFFFF;
	if((bank_len) > len){
		ClearErrorMessage();
		err_mess << "Segment length word exceeds valid words in buffer (" << bank_len << " > " << len << ")";
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	uint32_t type = (outbuff[0]>>16) & 0x0F;
	uint32_t Nwords = bank_len;
	uint32_t Nswapped = 1;
	switch(type){
		case 0x0a:  // 64 bit unsigned int
		case 0x08:  // 64 bit double
		case 0x09:  // 64 bit signed int
			swap_block((uint64_t*)&inbuff[1], Nwords/2, (uint64_t*)&outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x01:  // 32 bit unsigned int
		case 0x02:  // 32 bit float
		case 0x0b:  // 32 bit signed int
			swap_block(&inbuff[1], Nwords, &outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x05:  // 16 bit unsigned int
		case 0x04:  // 16 bit signed int
			swap_block((uint16_t*)&inbuff[1], Nwords*2, (uint16_t*)&outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x00:  // 32 bit unknown (not swapped)
		case 0x07:  // 8 bit unsigned int
		case 0x06:  // 8 bit signed int
			memcpy((uint8_t*)&outbuff[1], (uint8_t*)&inbuff[1], Nwords*sizeof(uint32_t));
			Nswapped += Nwords;
			break;
	}

	return Nswapped;
}

//---------------------------------
// swap_segment
//---------------------------------
uint32_t HDEVIO::swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO segment. 

	if(len < 1){
		SetErrorMessage("Attempt to swap segment with len<1");
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	// Swap header/length word
	swap_block(inbuff, 1, outbuff);
	uint32_t bank_len = outbuff[0] & 0xFFFF;
	if((bank_len) > len){
		ClearErrorMessage();
		err_mess << "Segment length word exceeds valid words in buffer (" << bank_len << " > " << len << ")";
		err_code = HDEVIO_BANK_TRUNCATED;
		Nerrors++;
		Nbad_events++;
		return 0;
	}
	
	uint32_t type = (outbuff[0]>>16) & 0x3F;
	uint32_t Nwords = bank_len;
	uint32_t Nswapped = 1;
	switch(type){
		case 0x0a:  // 64 bit unsigned int
		case 0x08:  // 64 bit double
		case 0x09:  // 64 bit signed int
			swap_block((uint64_t*)&inbuff[1], Nwords/2, (uint64_t*)&outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x01:  // 32 bit unsigned int
		case 0x02:  // 32 bit float
		case 0x0b:  // 32 bit signed int
			swap_block(&inbuff[1], Nwords, &outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x05:  // 16 bit unsigned int
		case 0x04:  // 16 bit signed int
			swap_block((uint16_t*)&inbuff[1], Nwords*2, (uint16_t*)&outbuff[1]);
			Nswapped += Nwords;
			break;
		case 0x00:  // 32 bit unknown (not swapped)
		case 0x07:  // 8 bit unsigned int
		case 0x06:  // 8 bit signed int
			memcpy((uint8_t*)&outbuff[1], (uint8_t*)&inbuff[1], Nwords*sizeof(uint32_t));
			Nswapped += Nwords;
			break;
	}

	return Nswapped;
}

//------------------------
// Print_fbuff
//------------------------
void HDEVIO::Print_fbuff(void)
{
	cout << endl;
	cout << "      fbuff: " << hex << (uint64_t)fbuff << dec << endl;
	cout << "      fnext: " << hex << (uint64_t)fnext << dec << endl;
	cout << "  fbuff_end: " << hex << (uint64_t)fbuff_end << dec << endl;
	cout << " fbuff_size: " << fbuff_size << endl;
	cout << "  fbuff_len: " << fbuff_len << endl;
	cout << "    _gcount: " << _gcount << endl;
}

//------------------------
// PrintEVIOBlockHeader
//------------------------
void HDEVIO::PrintEVIOBlockHeader(void)
{

	cout << endl;
	cout << "EVIO Block Header:" << endl;
	cout << "------------------------" << endl;
	cout << " Block Length: " << HexStr(buff[0]) << " (" << buff[0] << " words = " << (buff[0]>>(10-2)) << " kB)" << endl;
	cout << " Block Number: " << HexStr(buff[1]) << endl;
	cout << "Header Length: " << HexStr(buff[2]) << " (should be 8)" << endl;
	cout << "  Event Count: " << HexStr(buff[3]) << endl;
	cout << "   Reserved 1: " << HexStr(buff[4]) << endl;
	cout << "     Bit Info: " << HexStr(buff[5]>>8) << endl;
	cout << "      Version: " << HexStr(buff[5]&0xFF) << endl;
	cout << "   Reserved 3: " << HexStr(buff[6]) << endl;
	cout << "   Magic word: " << HexStr(buff[7]) << (swap_needed ? " (after swapping)":"") << endl;
	cout << "Byte swapping is" << (swap_needed ? " ":" not ") << "needed" << endl;
}

//------------------------
// PrintStats
//------------------------
void HDEVIO::PrintStats(void)
{
	cout << endl;
	cout << "EVIO Statistics for " << filename << " :" << endl;
	cout << "------------------------" << endl;
	cout << "    Nblocks: " << Nblocks << endl;
	cout << "    Nevents: " << Nevents << endl;
	cout << "    Nerrors: " << Nerrors << endl;
	cout << "Nbad_blocks: " << Nbad_blocks << endl;
	cout << "Nbad_events: " << Nbad_events << endl;
	cout << endl;
}


