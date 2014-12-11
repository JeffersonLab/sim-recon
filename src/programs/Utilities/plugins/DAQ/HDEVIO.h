// $Id$
//
//    File: HDEVIO.h
// Created: Wed Dec 10 07:22:00 EST 2014
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#ifndef _HDEVIO_
#define _HDEVIO_

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif


#include <stdint.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

class HDEVIO{
	public:
		HDEVIO(string filename);
		virtual ~HDEVIO();
		
		typedef struct{
			uint32_t length;
			uint32_t blocknum;
			uint32_t headerlen;
			uint32_t eventcnt;
			uint32_t reserved1;
			uint32_t bitinfo;
			uint32_t version;
			uint32_t reserved2;
			uint32_t magic;
		}BLOCKHEADER_t;
		
		enum{
			HDEVIO_OK=0,
			HDEVIO_EOF,
			HDEVIO_FILE_TRUNCATED,
			HDEVIO_EVENT_BIGGER_THAN_BLOCK,
			HDEVIO_BAD_BLOCK_HEADER,
			HDEVIO_BLOCKSIZE_GREATER_THAN_LIMIT,
			HDEVIO_MEMORY_ALLOCATION_ERROR,
			HDEVIO_NO_EVENTS_IN_BUFFER,
			HDEVIO_USER_BUFFER_TOO_SMALL,
			HDEVIO_BANK_TRUNCATED,
			HDEVIO_UNKNOWN_BANK_TYPE
		}ERRORCODE_t;
		
		string filename;
		ifstream ifs;
		bool is_open;
		
		uint32_t *buff;         // buffer holding current block (if any)
		uint32_t *next;         // Pointer to start of next EVIO event within buff
		uint32_t *buff_end;     // Pointer to word just past end of valid buffer words
		uint32_t buff_size;     // memory currently allocated for buff in words
		uint32_t buff_len;      // valid words in buff
		uint32_t buff_limit;    // maximum allowed allocation for buff_max
		bool swap_needed;       // true if block header indicates swapping is needed
		BLOCKHEADER_t *bh;      // =buff, but cast as a BLOCKHEADER_t*
		uint32_t last_event_len;// used to hold last event length in words if user buffer was
		                        // too small, this is how big is should be allocated
		
		stringstream err_mess;  // last error message
		uint32_t err_code;    // last error code
		
		uint64_t Nblocks;
		uint64_t Nevents;
		uint64_t Nerrors;
		uint64_t Nbad_blocks;
		uint64_t Nbad_events;

		bool ReadBlock(void);
		bool read(uint32_t *user_buff, uint32_t user_buff_len);

		uint32_t swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		uint32_t swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		uint32_t swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		void PrintEVIOBlockHeader(void);
		void PrintStats(void);

	protected:

		void ClearErrorMessage(void){ err_mess.str(""); err_mess.clear();}
		void SetErrorMessage(string mess){ ClearErrorMessage(); err_mess<<mess;}
	
	public:

		//------------------------
		// HexStr
		//------------------------
		inline string HexStr(uint32_t v)
		{
			char str[256];
			sprintf(str, "0x%08x", v);
			
			return string(str);
		}

		//------------------------
		// HexStr
		//------------------------
		inline string HexStr(uint16_t v)
		{
			char str[256];
			sprintf(str, "0x%04x", v);
			
			return string(str);
		}

		//------------------------
		// HexStr
		//------------------------
		inline string HexStr(uint8_t v)
		{
			char str[256];
			sprintf(str, "0x%02x", v);
			
			return string(str);
		}

		//---------------------------------
		// swap_block
		//---------------------------------
		inline void swap_block(uint16_t *inbuff, uint16_t len, uint16_t *outbuff)
		{
			for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
				uint16_t inword = *inbuff;  // copy word to allow using same buffer for input and output
				uint8_t *inptr  = (uint8_t*)&inword;
				uint8_t *outptr = (uint8_t*)&outbuff[1];
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
			}
		}

		//---------------------------------
		// swap_block
		//---------------------------------
		inline void swap_block(uint32_t *inbuff, uint32_t len, uint32_t *outbuff)
		{
			for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
				uint32_t inword = *inbuff;  // copy word to allow using same buffer for input and output
				uint8_t *inptr  = (uint8_t*)&inword;
				uint8_t *outptr = (uint8_t*)&outbuff[1];
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
			}
		}

		//---------------------------------
		// swap_block
		//---------------------------------
		inline void swap_block(uint64_t *inbuff, uint64_t len, uint64_t *outbuff)
		{
			for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
				uint64_t inword = *inbuff;  // copy word to allow using same buffer for input and output
				uint8_t *inptr  = (uint8_t*)&inword;
				uint8_t *outptr = (uint8_t*)&outbuff[1];
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;

				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
				*(--outptr) = *inptr++;
			}
		}

};

#endif // _HDEVIO_

