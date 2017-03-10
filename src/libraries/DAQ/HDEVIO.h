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
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

// ----- Stolen from evio.h -----------
#define swap64(x) ( (((x) >> 56) & 0x00000000000000FFL) | \
                         (((x) >> 40) & 0x000000000000FF00L) | \
                         (((x) >> 24) & 0x0000000000FF0000L) | \
                         (((x) >> 8)  & 0x00000000FF000000L) | \
                         (((x) << 8)  & 0x000000FF00000000L) | \
                         (((x) << 24) & 0x0000FF0000000000L) | \
                         (((x) << 40) & 0x00FF000000000000L) | \
                         (((x) << 56) & 0xFF00000000000000L) )

#define swap32(x) ( (((x) >> 24) & 0x000000FF) | \
                         (((x) >> 8)  & 0x0000FF00) | \
                         (((x) << 8)  & 0x00FF0000) | \
                         (((x) << 24) & 0xFF000000) )

#define swap16(x) ( (((x) >> 8) & 0x00FF) | \
                         (((x) << 8) & 0xFF00) )
//---------------------------------------

class HDEVIO{
	public:
		HDEVIO(string filename, bool read_map_file=true, int verbose=1);
		virtual ~HDEVIO();
		
		enum{
			HDEVIO_OK=0,
			HDEVIO_EOF,
			HDEVIO_FILE_NOT_OPEN,
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
		
		enum BLOCKTYPE{
			kBT_UNKNOWN,
			kBT_BOR,
			kBT_EPICS,
			kBT_SYNC,
			kBT_PRESTART,
			kBT_GO,
			kBT_PAUSE,
			kBT_END,
			kBT_PHYSICS
		};
		
		// The following structs are used to overlay the
		// the EVIO header and first few words. This is
		// used for scanning/mapping the file.
		typedef struct{
			uint32_t ts_header;
			uint32_t timestamp;
		}EPICSHEADER_t;

		typedef struct{
			uint32_t first_crate_len;
			uint32_t first_crate_header;
		}BORHEADER_t;

		typedef struct{
			uint32_t trigger_bank_len;
			uint32_t trigger_bank_header;
			uint32_t trigger_bank_segment_header;
			uint32_t first_event_hi; // n.b. contradicts documentation!
			uint32_t first_event_lo;
		}PHYSICSHEADER_t;
		
		typedef struct{
			// Standard EVIO 8 word block header
			uint32_t length;
			uint32_t blocknum;
			uint32_t headerlen;
			uint32_t eventcnt;
			uint32_t reserved1;
			uint32_t bitinfo;
			uint32_t reserved2;
			uint32_t magic;

			// first EVIO event in block
			uint32_t event_len;
			uint32_t header;

			// Next few words depend on type of data in block
			union {
				EPICSHEADER_t    epics;
				BORHEADER_t      bor;
				PHYSICSHEADER_t  physics;
			};
		}BLOCKHEADER_t;

		typedef struct{

			uint32_t event_len;
			uint32_t header;

			// Next few words depend on type of data in event
			union {
				EPICSHEADER_t    epics;
				BORHEADER_t      bor;
				PHYSICSHEADER_t  physics;
			};
		}EVENTHEADER_t;

		class EVIOEventRecord{
			public:
				streampos pos;
				uint32_t event_len;
				uint32_t event_header;
				uint64_t first_event;
				uint64_t last_event;
				BLOCKTYPE event_type;
		};

		class EVIOBlockRecord{
			public:
				streampos pos;
				uint32_t block_len;
				bool swap_needed;
				vector<EVIOEventRecord> evio_events;
				uint64_t first_event;
				uint64_t last_event;
				BLOCKTYPE block_type;
		};

		string filename;
		ifstream ifs;
		bool is_open;
		uint32_t *fbuff;
		uint32_t *fnext;
		uint32_t *fbuff_end;
		uint64_t fbuff_size;
		uint64_t fbuff_len;
		uint64_t _gcount;
		
		
		uint32_t *buff;           // buffer holding current block (if any)
		uint32_t *next;           // Pointer to start of next EVIO event within buff
		uint32_t *buff_end;       // Pointer to word just past end of valid buffer words
		uint32_t buff_size;       // memory currently allocated for buff in words
		uint32_t buff_len;        // valid words in buff
		uint32_t buff_limit;      // maximum allowed allocation for buff_max
		bool swap_needed;         // true if block header indicates swapping is needed
		BLOCKHEADER_t *bh;        // =buff, but cast as a BLOCKHEADER_t*
		streampos last_event_pos; // used to hold file position at last event read in
		uint32_t last_event_len;  // used to hold last event length in words if user buffer was
		                          // too small, this is how big is should be allocated
		
		int  VERBOSE;
		bool IGNORE_EMPTY_BOR;
		bool SKIP_EVENT_MAPPING;
		
		stringstream err_mess;  // last error message
		uint32_t err_code;    // last error code
		
		uint64_t Nblocks;
		uint64_t Nevents;
		uint64_t Nerrors;
		uint64_t Nbad_blocks;
		uint64_t Nbad_events;

		bool ReadBlock(void);
		bool read(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap=true);
		bool readSparse(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap=true);
		bool readNoFileBuff(uint32_t *user_buff, uint32_t user_buff_len, bool allow_swap=true);
		void rewind(void);

		uint32_t swap_bank(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		uint32_t swap_tagsegment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		uint32_t swap_segment(uint32_t *outbuff, uint32_t *inbuff, uint32_t len);
		void Print_fbuff(void);
		void PrintEVIOBlockHeader(void);
		void PrintStats(void);
		void PrintFileSummary(void);
		void SaveFileMap(string fname="");
		void ReadFileMap(string fname="", bool warn_if_not_found=false);

		uint32_t GetEventMask(void) { return event_type_mask; }
		uint32_t SetEventMask(uint32_t mask);
		uint32_t SetEventMask(string types_str);
		uint32_t AddToEventMask(string type_str);
		vector<EVIOBlockRecord>& GetEVIOBlockRecords(void);
		
	protected:
	
		uint32_t event_type_mask; 

		bool is_mapped;
		uint64_t total_size_bytes;
		vector<EVIOBlockRecord> evio_blocks;
		void MapBlocks(bool print_ticker=true);
		void MapEvents(BLOCKHEADER_t &bh, EVIOBlockRecord &br);
		vector<EVIOBlockRecord>::iterator sparse_block_iter;
		uint32_t sparse_event_idx;
		EVIOBlockRecord NB_block_record;
		streampos NB_next_pos;

		void ClearErrorMessage(void){ err_mess.str(""); err_mess.clear();}
		void SetErrorMessage(string mess){ ClearErrorMessage(); err_mess<<mess;}
		
		void buff_read(char* s, streamsize n);
		void buff_seekg (streamoff off, ios_base::seekdir way);
		streamsize buff_gcount() const { return (streamsize)_gcount; }
		
	
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
				*outbuff = swap16(*inbuff);
			}
		}

		//---------------------------------
		// swap_block
		//---------------------------------
		inline void swap_block(uint32_t *inbuff, uint32_t len, uint32_t *outbuff)
		{
			for(uint32_t i=0; i<len; i++, inbuff++, outbuff++){
				*outbuff = swap32(*inbuff);
			}
		}

		//---------------------------------
		// swap_block
		//---------------------------------
		inline void swap_block(uint64_t *inbuff, uint64_t len, uint64_t *outbuff)
		{
			for(uint64_t i=0; i<len; i++, inbuff++, outbuff++){
				*outbuff = swap64(*inbuff);
			}
		}

};

#endif // _HDEVIO_

