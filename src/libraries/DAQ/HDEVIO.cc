// $Id$
//
//    File: HDEVIO.cc
// Created: Wed Dec 10 07:22:00 EST 2014
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <unistd.h>
#include <cinttypes>
using namespace std;

#include "HDEVIO.h"

//---------------------------------
// HDEVIO    (Constructor)
//---------------------------------
HDEVIO::HDEVIO(string filename, bool read_map_file, int verbose):filename(filename),VERBOSE(verbose)
{
	// These must be initialized in case we return early
	// so they aren't deleted in the destructor if they
	// were never allocated.
	fbuff = NULL;
	buff  = NULL;

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
	last_event_pos = 0;
	last_event_len = 0;
	err_code = HDEVIO_OK;

	Nblocks = 0;
	Nevents = 0;
	Nerrors = 0;
	Nbad_blocks = 0;
	Nbad_events = 0;
	
	event_type_mask = 0xFFFF; // default to accepting all types
	is_mapped = false;
	
	NB_next_pos = 0;
	
	if(read_map_file) ReadFileMap(); // check if a map file exists and read it if it does
	
	IGNORE_EMPTY_BOR   = false;
	SKIP_EVENT_MAPPING = false;
	
	ifs.seekg(0, ios_base::end);
	total_size_bytes = ifs.tellg();
	ifs.seekg(0, ios_base::beg);
	
	is_open = true;
}

//---------------------------------
// ~HDEVIO    (Destructor)
//---------------------------------
HDEVIO::~HDEVIO()
{
	if(ifs.is_open()) ifs.close();
	if(buff ) delete[] buff;
	if(fbuff) delete[] fbuff;
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
bool HDEVIO::read(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap)
{
	/// Read the next EVIO event into the user supplied buffer.
	/// Return true if successful and false otherwise. Details of
	/// the error will be in err_mess.
	
	// If only certain event types are requested then
	// defer to the sparse reader
	if(event_type_mask != 0xFFFF){
		return readSparse(user_buff, user_buff_len, allow_swap);
	}

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
	if(swap_needed && allow_swap){
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
// readSparse
//---------------------------------
bool HDEVIO::readSparse(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap)
{
	/// This is an alternative to the read(...) method above that
	/// is used when the user has specified that only certain
	/// event types are desired. This really only makes sense
	/// for EPICS and SYNC events. This method does not use the
	/// fbuff system and instead reads directly from the file 
	/// after seeking to the desired location determined from
	/// a previously generated map.
	
	err_code = HDEVIO_OK;
	ClearErrorMessage();

	// Make sure we've mapped this file
	if(!is_mapped) MapBlocks();
	
	// Loop over all events of all blocks looking for the next
	// event matching the currently set type mask. 
	for(; sparse_block_iter!=evio_blocks.end(); sparse_block_iter++, sparse_event_idx = 0){

		// Filter out blocks of the wrong type
		EVIOBlockRecord &br = *sparse_block_iter;
		//		uint32_t type = (1 << br.block_type);

		for(; sparse_event_idx < br.evio_events.size(); sparse_event_idx++){
			EVIOEventRecord &er = sparse_block_iter->evio_events[sparse_event_idx];

			uint32_t etype = (1 << er.event_type);
			if( etype & event_type_mask ) break;
		}
		if(sparse_event_idx >= br.evio_events.size()) continue;
		
		EVIOEventRecord &er = sparse_block_iter->evio_events[sparse_event_idx];

		uint32_t event_len = er.event_len;
		last_event_len = event_len;
	
		// Check if user buffer is big enough to hold block
		if( event_len > user_buff_len ){
			ClearErrorMessage();
			err_mess << "user buffer too small for event (" << user_buff_len << " < " << event_len << ")";
			err_code = HDEVIO_USER_BUFFER_TOO_SMALL;
			return false;
		}

		// At this point we're committed to reading this event so go
		// ahead and increment pointer to next event so no matter
		// what happens below, we don't try reading it again.
		sparse_event_idx++;

		// Set file pointer to start of EVIO event (NOT block header!)
		last_event_pos = er.pos;
		ifs.seekg(last_event_pos, ios_base::beg);
		
		// Read data directly into user buffer
		ifs.read((char*)user_buff, event_len*sizeof(uint32_t));

		// Swap entire bank if needed
		swap_needed = br.swap_needed; // set flag in HDEVIO
		bool isgood = true;
		if(br.swap_needed && allow_swap){
			uint32_t Nswapped = swap_bank(user_buff, user_buff, event_len);
			isgood = (Nswapped == event_len);
		}
		
		// Double check that event length matches EVIO block header
		// but only if we either don't need to swap or need to and
		// were allowed to (otherwise, the test will almost certainly
		// fail!)
		if( (!br.swap_needed) || (br.swap_needed && allow_swap) ){
			if( (user_buff[0]+1) != event_len ){
				ClearErrorMessage();
				err_mess << "WARNING: EVIO bank indicates a different size than block header (" << event_len << " != " << (user_buff[0]+1) << ")";
				err_code = HDEVIO_EVENT_BIGGER_THAN_BLOCK;
				Nerrors++;
				Nbad_blocks++;
				return false;
			}
		}

		if(isgood) Nevents++;

		return isgood;
	}

	// If we got here then we did not find an event of interest
	// above. Report that there are no more events in the file.
	SetErrorMessage("No more events");
	err_code = HDEVIO_EOF;
	return false; // isgood=false
}

//---------------------------------
// readNoFileBuff
//---------------------------------
bool HDEVIO::readNoFileBuff(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap)
{
	/// This is an alternative to the read(...) method above that
	/// does not use a large primary file buffer. A single EVIO
	/// header is read in at a time and the events within the block
	/// mapped just like when using readSparse. The difference is
	/// that here, only a single block is mapped at a time rather
	/// than trying to map the entire file before starting. This
	/// gives a faster start up. This may be quicker for most
	/// desktop filesystems but may be slower for Lustre file systems
	/// that are configured for large volume data transfers and show
	/// perfomance degredation with small reads.
	
	err_code = HDEVIO_OK;
	ClearErrorMessage();
	
	// Check if we need to read in a block header using the 
	// current file position
	EVIOBlockRecord &br = NB_block_record;
	if(br.evio_events.empty()){

		// Check if we are at end of file
		uint64_t words_left_in_file = (total_size_bytes-NB_next_pos)/4;
		if( words_left_in_file == 8 ){ // (if <8 then let read below fail and return HDEVIO_FILE_TRUNCATED)
			SetErrorMessage("No more events");
			err_code = HDEVIO_EOF;
			return false;			
		}

		// read EVIO block header
		BLOCKHEADER_t bh;
		ifs.seekg( NB_next_pos, ios_base::beg);
		ifs.clear();
		ifs.read((char*)&bh, sizeof(bh));
		if(!ifs.good()){
			err_mess << "Error reading EVIO block header (truncated?)"<<endl;
			err_mess << "words_left_in_file (before read): " << words_left_in_file << endl;
			err_mess << "total_size_bytes: " << total_size_bytes << "   tellg: " << ifs.tellg();
			err_code = HDEVIO_FILE_TRUNCATED;
			return false;
		}

		// Check if we need to byte swap and simultaneously
		// verify header is good by checking magic word
		bool swap_needed = false;
		if(bh.magic==0x0001dac0){
			swap_needed = true;
		}else{
			if(bh.magic!=0xc0da0100){
				err_mess.str("Bad magic word");
				err_code = HDEVIO_BAD_BLOCK_HEADER;
				return false;
			}
		}
		
		if(swap_needed)swap_block((uint32_t*)&bh, sizeof(bh)>>2, (uint32_t*)&bh);
		
		Nblocks++;
		streampos pos = ifs.tellg() - (streampos)sizeof(bh);

		if( (uint64_t)(pos+(streampos)bh.length) >  total_size_bytes ){
			err_mess << "EVIO block extends past end of file!";
			err_code = HDEVIO_FILE_TRUNCATED;
			return false;
		}
		
		br.pos = pos;
		br.block_len = bh.length;
		br.swap_needed = swap_needed;
		br.first_event = 0;
		br.last_event = 0;
		
		MapEvents(bh, br);
		
		NB_next_pos = pos + (streampos)(bh.length<<2);
	}
	
	// Check if we did not find an event of interest above. 
	// If not, report that there are no more events in the file.
	if(br.evio_events.empty() || !ifs.good()){
		SetErrorMessage("No more events");
		err_code = HDEVIO_EOF;
		return false; // isgood=false
	}

	// Grab next event record
	EVIOEventRecord &er = br.evio_events.front();
	
	uint32_t event_len = er.event_len;
	last_event_len = event_len;

	// Check if user buffer is big enough to hold block
	if( event_len > user_buff_len ){
		ClearErrorMessage();
		err_mess << "user buffer too small for event (" << user_buff_len << " < " << event_len << ")";
		err_code = HDEVIO_USER_BUFFER_TOO_SMALL;
		return false;
	}

	// Set file pointer to start of EVIO event (NOT block header!)
	last_event_pos = er.pos;
	ifs.seekg(last_event_pos, ios_base::beg);
	
	// Read data directly into user buffer
	ifs.read((char*)user_buff, event_len*sizeof(uint32_t));
	if(!ifs.good()){
		SetErrorMessage("No more events");
		err_code = HDEVIO_EOF;
		return false; // isgood=false
	}

	// Remove EVIO Event record, effectively advancing to
	// next event for the next time we're called
	br.evio_events.erase(NB_block_record.evio_events.begin());

	// Swap entire bank if needed
	swap_needed = br.swap_needed; // set flag in HDEVIO
	bool isgood = true;
	if(br.swap_needed && allow_swap){
		uint32_t Nswapped = swap_bank(user_buff, user_buff, event_len);
		isgood = (Nswapped == event_len);
	}
	
	// Double check that event length matches EVIO block header
	// but only if we either don't need to swap or need to and
	// were allowed to (otherwise, the test will almost certainly
	// fail!)
	if( (!br.swap_needed) || (br.swap_needed && allow_swap) ){
		if( (user_buff[0]+1) != event_len ){
			ClearErrorMessage();
			err_mess << "WARNING: EVIO bank indicates a different size than block header (" << event_len << " != " << (user_buff[0]+1) << ")";
			err_code = HDEVIO_EVENT_BIGGER_THAN_BLOCK;
			Nerrors++;
			Nbad_blocks++;
			return false;
		}
	}

	if(isgood) Nevents++;

	return isgood;
}

//------------------------
// rewind
//------------------------
void HDEVIO::rewind(void)
{
	/// This can be used whe reading from a file to
	/// reset the file pointer and other position holders
	/// to the begining of the file. This is done when
	/// the "LOOP_FOREVER" option is used in the event source
	/// to continuously re-read a file, essesntially making
	/// it an infinite stream of events that can be used for
	/// testing.

	ifs.seekg(0, ios_base::beg);
	ifs.clear();
	
	sparse_block_iter = evio_blocks.begin();
	sparse_event_idx  = 0;
	
	NB_block_record.evio_events.clear();
	NB_next_pos = 0;
	
	ClearErrorMessage();
	err_code = HDEVIO_OK;
}

//------------------------
// SetEventMask
//------------------------
uint32_t HDEVIO::SetEventMask(uint32_t mask)
{
	uint32_t prev_mask = event_type_mask;
	event_type_mask = mask;

	return prev_mask;
}

//------------------------
// SetEventMask
//------------------------
uint32_t HDEVIO::SetEventMask(string types_str)
{
	uint32_t prev_mask = event_type_mask;

	event_type_mask = 0;
	if(types_str.find("BOR"    ) != string::npos) event_type_mask |= (1<<kBT_BOR);
	if(types_str.find("EPICS"  ) != string::npos) event_type_mask |= (1<<kBT_EPICS);
	if(types_str.find("PHYSICS") != string::npos) event_type_mask |= (1<<kBT_PHYSICS);

	return prev_mask;
}

//------------------------
// AddToEventMask
//------------------------
uint32_t HDEVIO::AddToEventMask(string type_str)
{
	return 0;
}

//------------------------
// GetEVIOBlockRecords
//------------------------
vector<HDEVIO::EVIOBlockRecord>& HDEVIO::GetEVIOBlockRecords(void)
{
	if(!is_mapped) MapBlocks();
	
	return evio_blocks;
}

//------------------------
// MapBlocks
//------------------------
void HDEVIO::MapBlocks(bool print_ticker)
{
	if(!is_open){
		err_mess.str("File is not open");
		err_code = HDEVIO_FILE_NOT_OPEN;
		return;
	}
	
	// Remember current file pos so we can restore it.
	streampos start_pos = ifs.tellg();
	
	if(print_ticker) cout << "Mapping EVIO file ..." << endl;
	
	// Rewind to beginning of file and loop over all blocks
	ifs.seekg(0, ios_base::beg);
	BLOCKHEADER_t bh;
	uint64_t Nblocks = 0;
	while(ifs.good()){
		ifs.read((char*)&bh, sizeof(bh));
		if(!ifs.good()) break;
		
		// Check if we need to byte swap and simultaneously
		// verify header is good by checking magic word
		bool swap_needed = false;
		if(bh.magic==0x0001dac0){
			swap_needed = true;
		}else{
			if(bh.magic!=0xc0da0100){
				err_mess.str("Bad magic word");
				err_code = HDEVIO_BAD_BLOCK_HEADER;
				EVIOBlockRecord br;
				br.pos = ifs.tellg() - (streampos)sizeof(bh);
				br.block_type = kBT_UNKNOWN;
				evio_blocks.push_back(br);
				break;
			}
		}
		
		if(swap_needed)swap_block((uint32_t*)&bh, sizeof(bh)>>2, (uint32_t*)&bh);
		
		Nblocks++;
		streampos block_len_bytes = (bh.length<<2)-sizeof(bh); // <<2 is for *4
		streampos pos = ifs.tellg() - (streampos)sizeof(bh);
		
		EVIOBlockRecord br;
		br.pos = pos;
		br.block_len = bh.length;
		br.swap_needed = swap_needed;
		br.first_event = 0;
		br.last_event = 0;

		// Categorize this block		
		uint32_t tag = bh.header>>16;
		uint32_t M   = bh.header&0xFF;
		
		switch(tag){
			case 0xFFD0: br.block_type = kBT_SYNC;      break;
			case 0xFFD1: br.block_type = kBT_PRESTART;  break;
			case 0xFFD2: br.block_type = kBT_GO;        break;
			case 0xFFD3: br.block_type = kBT_PAUSE;     break;
			case 0xFFD4: br.block_type = kBT_END;       break;
			case 0x0060: br.block_type = kBT_EPICS;     break;
			case 0x0070: br.block_type = kBT_BOR;       break;
			case 0xFF50:
			case 0xFF51:
			case 0xFF70:
				br.block_type   = kBT_PHYSICS;
				br.first_event  = bh.physics.first_event_lo;
				br.first_event += ((uint64_t)bh.physics.first_event_hi)<<32;
				br.last_event   = br.first_event + (uint64_t)M - 1;
				break;
			default:
				br.block_type   = kBT_UNKNOWN;
				_DBG_ << "Uknown tag: " << hex << tag << dec << endl;
		}
		
		// Scan through and map all events within this block
		if( !SKIP_EVENT_MAPPING ) MapEvents(bh, br);

		// Add block to list
		evio_blocks.push_back(br);
		
		// Update ticker
		if(print_ticker){
			if((Nblocks%500) == 0){
				uint64_t total_MB = total_size_bytes>>20;
				uint64_t read_MB  = ifs.tellg()>>20;
				if(Nblocks==0) cout << endl;
				cout << Nblocks << " blocks scanned (" << read_MB << "/" << total_MB << " MB " << (100*read_MB/total_MB) << "%)     \r";
				cout.flush();
			}
		}

		// Advance file pointer to start of next EVIO block header
		ifs.seekg(block_len_bytes, ios_base::cur);
	}
	
	if(print_ticker) cout << endl;
	
	// Setup iterators for sparse reading
	sparse_block_iter = evio_blocks.begin();
	sparse_event_idx = 0;

	// Restore file pos and set flag that file has been mapped
	ifs.clear();
	ifs.seekg(start_pos, ios_base::beg);
	is_mapped = true;
}

//---------------------------------
// MapEvents
//---------------------------------
void HDEVIO::MapEvents(BLOCKHEADER_t &bh, EVIOBlockRecord &br)
{
	/// This is called if the EVIO  block header indicates that
	/// it contains more than one top-level event. The position
	/// and event length of the first top-level event are passed
	/// in as starting parameters.
	
	// Record stream position upon entry so we can restore it at end
	streampos start_pos = ifs.tellg();
	
	// Calculate stream position of EVIO event
	streampos pos = start_pos -(streampos)sizeof(BLOCKHEADER_t)  + (streampos)(8<<2); // (8<<2) is 8 word EVIO block header times 4bytes/word 
	ifs.seekg(pos, ios_base::beg);
	
	EVENTHEADER_t myeh;
	for(uint32_t i=0; i<bh.eventcnt; i++){
	
		// For the first iteration through this loop,
		// we use the event header that has already been
		// read in as part of the block header.
		EVENTHEADER_t *eh = (EVENTHEADER_t*)&bh.event_len;
		if(i!=0){
			// Read in first few words of event
			eh = &myeh;
			ifs.read((char*)eh, sizeof(EVENTHEADER_t));
			if(!ifs.good()) break;
			if(br.swap_needed)swap_block((uint32_t*)eh, sizeof(EVENTHEADER_t)>>2, (uint32_t*)eh);
		}else{
			ifs.seekg(sizeof(EVENTHEADER_t), ios_base::cur);
		}

		if (eh->event_len < 2) {
			// Before disabling this warning (or hiding it behind a VERBOSE flag)
			// you should ask yourself the question, "Is this something that we
			// should simply be ignoring, garbage bytes in the input evio file?"
			std::cout << "HDEVIO::MapEvents warning - " << "Attempt to swap bank with len<2 (len="<<eh->event_len<<" header="<<hex<<eh->header<<dec<<" pos=" << pos << " tellg=" << ifs.tellg() << " i=" << i << ")" << std::endl;
			
			// Reference run 20495: Seems ROL is putting BOR bank header, but
			// the bank has no data in it. For this case we can go forward, but
			// this is really a problem for production data. Detect this and warn
			// user.
			if(eh->event_len==1 && (eh->header&0xFFFF00FF)==0x00700001){
				_DBG__;
				_DBG_ << "WARNING: This looks like an empty BOR event. BOR configuration" << endl;
				_DBG_ << "         data will not be available and it is unlikely you will" << endl;
				_DBG_ << "         be able to do anything beyond the digihit level. " << endl;
				_DBG__;
				
				if(IGNORE_EMPTY_BOR){
					EVIOEventRecord er;
					er.pos = pos;
					er.event_len   = eh->event_len + 1; // +1 to include length word
					er.event_type  = kBT_UNKNOWN;
					er.first_event = 0;
					er.last_event  = 0;
					br.evio_events.push_back(er);
				}else{
					_DBG_ << "         The program will (probably) stop now." << endl;
					_DBG_ << "         To avoid stopping, re-run with EVIO:IGNORE_EMPTY_BOR=1 ." << endl;				}
			}
			
			Nbad_events++;
			Nerrors++;
			// --i;  // This caused an infinite loop when reading hd_rawdata_020058_000.evio DL
			streampos delta = (streampos)((eh->event_len+1)<<2) - 
                              (streampos)sizeof(EVENTHEADER_t);
			ifs.seekg(delta, ios_base::cur);
			pos += (streampos)((eh->event_len+1)<<2);
			continue;
		}

		EVIOEventRecord er;
		er.pos = pos;
		er.event_len   = eh->event_len + 1; // +1 to include length word
		er.event_header= eh->header;
		er.event_type  = kBT_UNKNOWN;
		er.first_event = 0;
		er.last_event  = 0;

		uint32_t tag = eh->header>>16;
		uint32_t M   = eh->header&0xFF;
		
		switch(tag){
			case 0xFFD0: er.event_type = kBT_SYNC;       break;
			case 0xFFD1: er.event_type = kBT_PRESTART;   break;
			case 0xFFD2: er.event_type = kBT_GO;         break;
			case 0xFFD3: er.event_type = kBT_PAUSE;      break;
			case 0xFFD4: er.event_type = kBT_END;        break;
			case 0x0060: er.event_type = kBT_EPICS;      break;
			case 0x0070: er.event_type = kBT_BOR;        break;
			case 0xFF50:
			case 0xFF51:
			case 0xFF70:
				er.event_type = kBT_PHYSICS;
				er.first_event  = eh->physics.first_event_lo;
				er.first_event += ((uint64_t)eh->physics.first_event_hi)<<32;
				er.last_event   = er.first_event + (uint64_t)M - 1;
				if(er.first_event < br.first_event) br.first_event = er.first_event;
				if(er.last_event  > br.last_event ) br.last_event  = er.last_event;
				break;
			default:
				if(VERBOSE>1) _DBG_ << "Uknown tag: " << hex << tag << dec << endl;
		}
		
		br.evio_events.push_back(er);
		
		// Move file position to start of next event
		streampos delta = (streampos)((eh->event_len+1)<<2) - (streampos)sizeof(EVENTHEADER_t);
		ifs.seekg(delta, ios_base::cur);
		pos += (streampos)((eh->event_len+1)<<2);
	}

	ifs.clear();
	ifs.seekg(start_pos, ios_base::beg);
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
			if( inbuff!=outbuff ) memcpy((uint8_t*)&outbuff[2], (uint8_t*)&inbuff[2], Nwords*sizeof(uint32_t));
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
			if( inbuff!=outbuff ) memcpy((uint8_t*)&outbuff[1], (uint8_t*)&inbuff[1], Nwords*sizeof(uint32_t));
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
	uint64_t Nblocks = this->Nblocks;
	uint64_t Nevents = this->Nevents;
	
	if(is_mapped){
		Nblocks = evio_blocks.size();
		Nevents = 0;
		for(auto b : evio_blocks) Nevents += b.evio_events.size();
	}

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

//------------------------
// PrintFileSummary
//------------------------
void HDEVIO::PrintFileSummary(void)
{
	if(!is_mapped) MapBlocks();

	uint32_t Nsync     = 0;
	uint32_t Nprestart = 0;
	uint32_t Ngo       = 0;
	uint32_t Npause    = 0;
	uint32_t Nend      = 0;
	uint32_t Nepics    = 0;
	uint32_t Nbor      = 0;
	uint32_t Nphysics  = 0;
	uint32_t Nunknown  = 0;
	uint32_t Nblockunknown = 0;

	uint64_t first_event = 0;
	uint64_t last_event  = 0;

	uint32_t map_size = evio_blocks.size()*sizeof(EVIOBlockRecord);
	
	set<uint32_t> block_levels;
	set<uint32_t> events_in_block;

	// Loop over all EVIO block records
	for(uint32_t i=0; i<evio_blocks.size(); i++){
		EVIOBlockRecord &br = evio_blocks[i];
		events_in_block.insert(br.evio_events.size());
		map_size += br.evio_events.size()*sizeof(EVIOEventRecord);

		uint32_t Nunknown_prev = Nunknown;
		for(uint32_t j=0; j<br.evio_events.size(); j++){
			EVIOEventRecord &er = br.evio_events[j];
			
			uint32_t block_level;
			switch(er.event_type){
				case kBT_UNKNOWN:    Nunknown++;     break;
				case kBT_SYNC:       Nsync++;        break;
				case kBT_PRESTART:   Nprestart++;    break;
				case kBT_GO:         Ngo++;          break;
				case kBT_PAUSE:      Npause++;       break;
				case kBT_END:        Nend++;         break;
				case kBT_EPICS:      Nepics++;       break;
				case kBT_BOR:        Nbor++;         break;
				case kBT_PHYSICS:
					block_level = (uint32_t)((er.last_event - er.first_event) + 1);
					block_levels.insert(block_level);
					Nphysics += block_level;
					if(er.first_event<first_event || first_event==0) first_event = er.first_event;
					if(er.last_event>last_event) last_event = er.last_event;
					break;
				default:
					break;
			}
			
			//_DBG_ << "Block " << i << "  event " << j << "  " << er.last_event <<" - " << er.first_event << " = " << block_level << endl;
		}
		
		if( (Nunknown-Nunknown_prev) > 0 ) Nblockunknown++;
	}
	
	// form succint string of block levels
	stringstream ss;
	set<uint32_t>::iterator it = block_levels.begin();
	for(; it!=block_levels.end(); it++) ss << *it << ",";
	string sblock_levels = ss.str();
	if(!sblock_levels.empty()) sblock_levels.erase(sblock_levels.length()-1);

	// form succint string of events per block
	ss.str("");
	it = events_in_block.begin();
	for(; it!=events_in_block.end(); it++){
		uint32_t val = *it;
		ss << val;
		if(++it==events_in_block.end()) break;
		if( *it == ++val ){
			ss << "-";
			for(it++; it!=events_in_block.end(); it++){
				if( *it != ++val ){
					ss << (val-1) << ",";
					it--;
					break;
				}
			}
		}else{
			ss << ",";
		}
		if( it==events_in_block.end() ) break;
	}
	string sevents_in_block = ss.str();
	
	// Print results
	PrintStats();
	cout << "EVIO file size: " << (total_size_bytes>>20) << " MB" <<endl;
	cout << "EVIO block map size: " << (map_size>>10) << " kB" <<endl;
	cout << "first event: " << first_event << endl;
	cout << "last event: " << last_event << endl;

	cout << endl;
	cout << "             block levels = " << sblock_levels << endl;
	cout << "         events per block = " << sevents_in_block << endl;
	cout << "                    Nsync = " << Nsync << endl;
	cout << "                Nprestart = " << Nprestart << endl;
	cout << "                      Ngo = " << Ngo << endl;
	cout << "                   Npause = " << Npause << endl;
	cout << "                     Nend = " << Nend << endl;
	cout << "                   Nepics = " << Nepics << endl;
	cout << "                     Nbor = " << Nbor << endl;
	cout << "                 Nphysics = " << Nphysics << endl;
	cout << "                 Nunknown = " << Nunknown << endl;
	cout << " blocks with unknown tags = " << Nblockunknown << endl;
	cout << endl;
}

//------------------------
// SaveFileMap
//------------------------
void HDEVIO::SaveFileMap(string fname)
{
	// Make sure file has been mapped
	if(!is_mapped) MapBlocks();
	
	// Open output file
	if(fname=="") fname = filename + ".map";
	ofstream ofs(fname.c_str());
	if(!ofs.is_open()){
		cerr << "Unable to open \""<<fname<<"\" for writing!" << endl;
		return;
	}
	
	cout << "Writing EVIO file map to: " << fname << endl;
	
	char str[256];
	time_t t = time(NULL);
	struct tm *tmp = localtime(&t);
	strftime(str, 255, "%c", tmp);
	
	ofs << "#" << endl;
	ofs << "# HDEVIO block map for " << filename << endl;
	ofs << "#" << endl;
	ofs << "# generated " << str << endl;
	ofs << "#" << endl;
	ofs << "# " << endl;
	ofs << "swap_needed: " << swap_needed << endl;
	ofs << "--------- Start of block data --------" << endl;
	ofs << "#  pos    block_len  first_evt last_evt block_type" << endl;
	ofs << "#  + pos    evt_len  evt_header first_evt last_evt event_type" << endl;

	for(auto br : evio_blocks){
		char line[512];
		sprintf(line, "0x%08x 0x%06x  %8" PRIu64 "  %8" PRIu64 "     %d", (uint32_t)br.pos, br.block_len, br.first_event, br.last_event, br.block_type);
		ofs << line << endl;
		
		for(auto er : br.evio_events){
			sprintf(line, "+ 0x%08x 0x%06x 0x%08x %8" PRIu64 " %8" PRIu64 " %d", (uint32_t)er.pos, er.event_len, er.event_header, er.first_event, er.last_event, er.event_type);
			ofs << line << endl;
		}
	}
	ofs << "# --- End of map ---" << endl;
	// Close output file	
	ofs.close();
	
	cout << "Done" << endl;
}

//------------------------
// ReadFileMap
//------------------------
void HDEVIO::ReadFileMap(string fname, bool warn_if_not_found)
{
	// Open input file
	if(VERBOSE>4) cout << " Attempting to read EVIO map file \"" << fname << "\" for \"" << filename << "\"" << endl;
	if(fname=="") {
		
		// No map file name given. Form a list of potential ones
		// in the order they should be checked.
		string dname = ".";
		string bname = filename;
		auto pos = filename.find_last_of("/");
		if(pos != string::npos){
			dname = filename.substr(0, pos);
			bname = filename.substr(pos+1, filename.size()-pos);
		}

		vector<string> fnames;
		fnames.push_back(filename + ".map");
		fnames.push_back(dname + "/filemaps/" + bname + ".map");
		fnames.push_back(bname + ".map");
		fnames.push_back(filename + ".bmap");
		fnames.push_back(dname + "/filemaps/" + bname + ".bmap");
		fnames.push_back(bname + ".bmap");
		
		// Loop over possible names until we find one that is readable
		for(string f : fnames){
			if(VERBOSE>2) cout << "Checking for EVIO map file: " << f << " ...";
			if( access(f.c_str(), R_OK) != 0 ) {
				if(VERBOSE>2) cout << "no" << endl;
				continue;
			}
			if(VERBOSE>2)cout << "yes" << endl;
			fname = f;
			break;
		}		
	}
	
	if(fname=="") return;
	
	// Open map file
	ifstream ifs(fname.c_str());
	if(!ifs.is_open()){
		if(warn_if_not_found) cerr << "Unable to open \""<<fname<<"\" for reading!" << endl;
		return;
	}

	// Check if file was closed cleanly
	string eof_string("# --- End of map ---");
	ifs.seekg(-eof_string.length()*2, ifs.end); // not guaranteed how many char's endl is
	bool closed_cleanly = false;
	string line;
	while(getline(ifs,line)){
		if(string(line).find(eof_string)!=string::npos) closed_cleanly = true;
	}
	if(!closed_cleanly){
		cerr << "Found map file \"" << fname << "\" but it wasn't closed cleanly. Ignoring." << endl;
		ifs.close();
		return;
	}
	
	// Reset file pointers to start of file
	ifs.clear();
	ifs.seekg(0);
	
	// Loop over header
	while(getline(ifs, line)){
		
		if(line.length() < 5   ) continue;
		if(line.find("#") == 0 ) continue;
		if(line.find("Start of block data") != string::npos) break;
	}
	
	// Loop over body
	EVIOBlockRecord br;
	bool first_block_found = false;
	string s;
	while(getline(ifs, s)){
		stringstream ss(s);
		if(ss.str().find(eof_string) != string::npos ) break; // end of map trailer
		if(ss.str().find("#") == 0 ) continue; // ignore all other comment lines
		
		if(ss.str().find("+") == 0 ){
			// EVIO Event Record
			string tmp;
			uint64_t tmp64;
			EVIOEventRecord er;
			ss >> tmp; // '+"
			ss << hex;
			ss >> tmp64; er.pos = tmp64; // operator>> won't stream directly to streampos
			ss >> er.event_len;
			ss >> er.event_header;
			ss << dec;
			ss >> er.first_event;
			ss >> er.last_event;
			ss >> tmp64; er.event_type = (BLOCKTYPE)tmp64; // operator>> won't stream directly to BLOCKTYPE
			br.evio_events.push_back(er);
		}else{
			// EVIO Block Record
			if(first_block_found){
				evio_blocks.push_back(br);
			}else{
				first_block_found = true;
			}
			br.evio_events.clear();
			uint64_t tmp64;
			ss << hex;
			ss >> tmp64; br.pos = tmp64; // operator>> won't stream directly to streampos
			ss >> br.block_len;
			ss << dec;
			ss >> br.first_event;
			ss >> br.last_event;
			ss >> tmp64; br.block_type = (BLOCKTYPE)tmp64; // operator>> won't stream directly to BLOCKTYPE
		}
	}
	if(!br.evio_events.empty()) evio_blocks.push_back(br);
	
	is_mapped = true;
	cout << "Read EVIO file map from: " << fname << endl;
}

