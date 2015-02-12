//
// $Id: $
// $HeadURL: $
//

#include <stdlib.h>
#include <stdint.h>

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#ifndef _DBG_
#define _DBG_  cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

#define EVIO_SWAP64(x) ( (((x) >> 56) & 0x00000000000000FFL) | \
                         (((x) >> 40) & 0x000000000000FF00L) | \
                         (((x) >> 24) & 0x0000000000FF0000L) | \
                         (((x) >> 8)  & 0x00000000FF000000L) | \
                         (((x) << 8)  & 0x000000FF00000000L) | \
                         (((x) << 24) & 0x0000FF0000000000L) | \
                         (((x) << 40) & 0x00FF000000000000L) | \
                         (((x) << 56) & 0xFF00000000000000L) )

#define EVIO_SWAP32(x) ( (((x) >> 24) & 0x000000FF) | \
                         (((x) >> 8)  & 0x0000FF00) | \
                         (((x) << 8)  & 0x00FF0000) | \
                         (((x) << 24) & 0xFF000000) )

#define EVIO_SWAP16(x) ( (((x) >> 8) & 0x00FF) | \
                         (((x) << 8) & 0xFF00) )

string FILENAME = "";
bool PRINT_BLOCK_HEADERS = false;
bool PRINT_BANK_INFOS = false;

bool swap_needed = false;

uint32_t Nblocks_read = 0;
uint32_t Nevents_read = 0;
uint32_t Nbanks_read = 0;
uint32_t Nerrors = 0;

string HexStr(uint32_t v);
string HexStr(uint16_t v);
string HexStr(uint8_t v);
void PrintEVIOBlockHeader(uint32_t *buff);
void PrintEVIOBankInfo(uint32_t *buff);
uint32_t ProcessBlock(uint32_t blknum, uint32_t *istart, uint32_t *iend, uint32_t &Ntoplevelbanks);
void ParseCommandLineArguments(int narg, char *argv[]);
void Usage(void);

//------------------------
// main
//------------------------
int main(int narg, char *argv[])
{
	ParseCommandLineArguments(narg, argv);
	
	cout << endl;

	// Open EVIO file
	ifstream ifs(FILENAME.c_str());
	if(!ifs.is_open()){
		cerr << "Unable to open file: " << FILENAME << " !!" << endl;
		exit(-1);
	}
	cout << "Opened file: " << FILENAME << endl;
	
	// Allocate initial buffer
	uint32_t buff_size_words = 1024;
	uint32_t *buff = new uint32_t[buff_size_words];

	while(ifs.good()){

	// Read in EVIO block header
		ifs.read((char*)buff, 8*sizeof(uint32_t));
		uint32_t valid_words = ifs.gcount()/sizeof(uint32_t);
		if(valid_words != 8){
			cerr << "Could not read in 8 word EVIO block header!" << endl;
			break;
			//exit(-2);
		}

		// Check endianess
		if(buff[7]!=0xc0da100 && buff[7]!=0x0001dac0){
			cerr << "Magic word not valid!: " << HexStr(buff[7]) << endl;
			exit(-3);
		}
		swap_needed = (buff[7]==0x0001dac0);

		// Swap header (if needed)
		if(swap_needed) for(uint32_t i=0; i<8; i++) buff[i] = EVIO_SWAP32(buff[i]);

		// Print block header
		if(PRINT_BLOCK_HEADERS) PrintEVIOBlockHeader(buff);

		// Re-allocate buffer if needed so we can read in entire block
		uint32_t block_length = buff[0];
		if(buff_size_words < block_length){
			if(buff) delete[] buff;
			buff_size_words = block_length;
			buff = new uint32_t[buff_size_words];
		}

		// Set file pointer back to start of EVIO block header and read entire block in
		ifs.seekg(-8*sizeof(uint32_t), ifs.cur);
		ifs.read((char*)buff, block_length*sizeof(uint32_t));
		valid_words = ifs.gcount()/sizeof(uint32_t);
		if(valid_words < block_length){
			cerr << "Error reading in EVIO block!";
			if(valid_words>=2) cerr << "  (block number: " << (swap_needed ? EVIO_SWAP32(buff[1]):buff[1]) << " )";
			cerr << endl;
			cerr << "block size is " << block_length << " words but only able to read " << valid_words << " words!" << endl;
			Nerrors++;
			break;
			//exit(-4);
		}
		if(swap_needed) for(uint32_t i=0; i<8; i++) buff[i] = EVIO_SWAP32(buff[i]); // Swap EVIO block header

		uint32_t Ntoplevelbanks;
		ProcessBlock(buff[1], &buff[8], &buff[valid_words], Ntoplevelbanks);
		
		Nblocks_read++;
		Nevents_read += buff[3];

		static time_t last_time = time(NULL);
		time_t now = time(NULL);
		if(now != last_time){
			last_time = now;
			cout << "read " << Nblocks_read << " EVIO blocks   \r";
			cout.flush();
		}
	
	} // while(ifs.good())
	cout << endl;
	
	// Close EVIO file
	ifs.close();
	
	if(buff) delete[] buff;

	cout << endl;
	cout << "EVIO blocks read: " << Nblocks_read << endl;
	cout << "    Events Found: " << Nevents_read << endl;
	cout << "          Errors: " << Nerrors << endl;

	return 0;
}

//------------------------
// HexStr
//------------------------
string HexStr(uint32_t v)
{
	char str[256];
	sprintf(str, "0x%08x", v);
	
	return string(str);
}

//------------------------
// HexStr
//------------------------
string HexStr(uint16_t v)
{
	char str[256];
	sprintf(str, "0x%04x", v);
	
	return string(str);
}

//------------------------
// HexStr
//------------------------
string HexStr(uint8_t v)
{
	char str[256];
	sprintf(str, "0x%02x", v);
	
	return string(str);
}

//------------------------
// PrintEVIOBlockHeader
//------------------------
void PrintEVIOBlockHeader(uint32_t *buff)
{
	uint32_t first_event = Nevents_read + 1;
	uint32_t last_event  = first_event + buff[3] - 1;

	cout << endl;
	cout << "EVIO Block Header: " << Nblocks_read << " (events " << first_event << " - " << last_event << ")" << endl;
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
// PrintEVIOBankInfo
//------------------------
void PrintEVIOBankInfo(uint32_t *buff)
{
	uint16_t tag  = buff[1]>>16;
	uint8_t num  = buff[1]&0xFF;
	uint8_t type = (buff[1]>>8)&0xFF;
	string typestr = "unknown";
	switch(type){
		case 0x00: typestr = "32-bit unknown (not swapped)";  break;
		case 0x01: typestr = "32-bit unsigned int";  break;
		case 0x02: typestr = "32-bit float";  break;
		case 0x03: typestr = "8-bit char*";  break;
		case 0x04: typestr = "16-bit signed short";  break;
		case 0x05: typestr = "16-bit unsigned short";  break;
		case 0x06: typestr = "8-bit signed char";  break;
		case 0x07: typestr = "8-bit unsigned char";  break;
		case 0x08: typestr = "64-bit double";  break;
		case 0x09: typestr = "64-bit signed int";  break;
		case 0x0a: typestr = "64-bit unsigned int";  break;
		case 0x0b: typestr = "32-bit signed int";  break;
		case 0x0c: typestr = "TAGSEGMENT";  break;
		case 0x0d: typestr = "SEGMENT";  break;
		case 0x0e: typestr = "BANK";  break;
		case 0x0f: typestr = "COMPOSITE";  break;
		case 0x10: typestr = "BANK";  break;
		case 0x20: typestr = "SEGMENT";  break;
		case 0x21: typestr = "Hollerit*";  break;
		case 0x22: typestr = "N value*";  break;
		default:   typestr = "unknown"; break;
	}

	cout << endl;
	cout << "EVIO Bank Info:" << endl;
	cout << "------------------------" << endl;
	cout << "Bank Length: " << HexStr(buff[0]) << " (" << buff[0] << " words = " << (buff[0]>>(10-2)) << " kB)" << endl;
	cout << "    Tag/Num: " << HexStr(buff[1]) << " tag=" << HexStr(tag) << " num=" << HexStr(num) << endl;
	cout << "       Type: "  << HexStr(type) << " (" << typestr << ")" << endl;

	cout << endl;
}

//------------------------
// ProcessBlock
//------------------------
uint32_t ProcessBlock(uint32_t blknum, uint32_t *istart, uint32_t *iend, uint32_t &Ntoplevelbanks)
{
	uint32_t *iptr = istart;

	uint32_t Nexpected_words = ((unsigned long)iend - (unsigned long)istart)/sizeof(uint32_t);
	
	Ntoplevelbanks = 0;
	uint32_t Ntot_words = 0;
	while(iptr < iend){
		if(swap_needed) for(uint32_t i=0; i<2; i++) iptr[i] = EVIO_SWAP32(iptr[i]);
		
		if(PRINT_BANK_INFOS) PrintEVIOBankInfo(iptr);
		
		uint32_t length = iptr[0];
		uint32_t header = iptr[1];

		iptr = &iptr[length +1];
		Ntoplevelbanks++;
		Ntot_words += length+1;
	}
	
	if(Nexpected_words != Ntot_words){
		static uint32_t Nwarnings=0;
		Nerrors++;
		Nwarnings++;
		if(Nwarnings<=10) _DBG_<< "WARNING! block " << blknum << " : Nexpected_words="<<Nexpected_words<<" while Ntot_words="<<Ntot_words<<endl;
		if(Nwarnings==10){
			cout << "At least 10 blocks seen where the word count obtained from" << endl;
			cout << "bank headers does not match the word count from the EVIO" << endl;
			cout << "block header! Now more warnings will be printed." << endl;
			//exit(-1);
		}
	}
	
	return iptr - istart;
}

//------------------------
// ParseCommandLineArguments
//------------------------
void ParseCommandLineArguments(int narg, char *argv[])
{
	for(int i=1; i<narg; i++){
		string arg = argv[i];
		string next = (i+1)<narg ? argv[i]:"";
		bool used_next = false;
		
		if(arg == "-h"    ) Usage();
		if(arg == "--help") Usage();
		if(arg[0] != '-') FILENAME = arg;
		
		if(used_next) i++;
	}
	
	if(FILENAME.length() == 0) Usage();
}

//------------------------
// Usage
//------------------------
void Usage(void)
{
	cout << endl;
	cout << "Usage:" << endl;
	cout << "     evio_check file.evio" << endl;
	cout << endl;
	cout << endl;
	
	exit(0);
}
