
#include <swap_bank.h>

#include <string.h>

#include <iostream>
using namespace std;

#include <JANA/JException.h>
using namespace jana;

//---------------------------------
// swap_bank
//---------------------------------
uint32_t swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO bank. If the bank contains data, it is automatically
	/// swapped according to it's type. If the bank is a container of other
	/// containers, then this repeatedly calls the swapper methods for the
	/// appropriate container type (bank, tagsegment, segment). This means
	/// that this method will be recursive in the cases where it is a bank
	/// of banks.

	if(len < 2) throw JException("Attempt to swap bank with len<2", __FILE__, __LINE__);
	
	// Swap length and header words
	swap_block(inbuff, 2, outbuff);
	uint32_t bank_len = outbuff[0];
	if((bank_len+1) > len){
		stringstream ss;
		ss << "WARNING: Bank length word exceeds valid words in buffer (" << bank_len+1 << " > " << len << ")";
		throw JException(ss.str(), __FILE__, __LINE__);
	}
	if( bank_len < 1 ){
		throw JException("EVIO bank length word is zero in swap_bank!", __FILE__, __LINE__);
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
			stringstream ss;
			ss << "WARNING: unknown bank type (0x" << hex << type << dec << ")";
			throw JException(ss.str(), __FILE__, __LINE__);

			break;
	}

	return Nswapped;
}

//---------------------------------
// swap_tagsegment
//---------------------------------
uint32_t swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO tagsegment. 

	if(len < 1) throw JException("Attempt to swap segment with len<1");
	
	// Swap header/length word
	swap_block(inbuff, 1, outbuff);
	uint32_t bank_len = outbuff[0] & 0xFFFF;
	if((bank_len) > len){
		stringstream ss;
		ss << "Segment length word exceeds valid words in buffer (" << bank_len << " > " << len << ")";
		throw JException(ss.str(), __FILE__, __LINE__);
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
			if( inbuff!=outbuff ) memcpy((uint8_t*)&outbuff[1], (uint8_t*)&inbuff[1], Nwords*sizeof(uint32_t));
			Nswapped += Nwords;
			break;
	}

	return Nswapped;
}

//---------------------------------
// swap_segment
//---------------------------------
uint32_t swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len)
{
	/// Swap an EVIO segment. 

	if(len < 1) throw JException("Attempt to swap segment with len<1");
	
	// Swap header/length word
	swap_block(inbuff, 1, outbuff);
	uint32_t bank_len = outbuff[0] & 0xFFFF;
	if((bank_len) > len){
		stringstream ss;
		ss << "Segment length word exceeds valid words in buffer (" << bank_len << " > " << len << ")";
		throw JException(ss.str(), __FILE__, __LINE__);
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
