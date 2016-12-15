// $Id$
//
//    File: HDET.h
// Created: Fri Apr 22 07:05:40 EDT 2016
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#ifndef _HDET_
#define _HDET_

#include <list>
#include <string>
#include <utility>
#include <sstream>
#include <iostream>
using namespace std;

#include <JANA/jerror.h>

#ifdef HAVE_ET
#include <et.h>
#endif // HAVE_ET


class HDET{
	public:

		enum{
			HDET_OK=0,            // Everything is hunky-dory
			HDET_ERROR,           // General error
			HDET_NO_ET_SUPPORT,   // ET library not compiled in
			HDET_TIMEOUT,         // Timeout trying to read ET event
			HDET_ET_SYSTEM_DEAD,  // ET system is dead
			HDET_ALLOC_FAILED,    // Failed to allocate a buffer
			HDET_BAD_FORMAT       // EVIO data in ET event has bad format
		}ERRORCODE_t;


		         HDET(string source_name, int ET_STATION_NEVENTS=10, bool ET_STATION_CREATE_BLOCKING=false);
		virtual ~HDET();

		    bool read(uint32_t* &buff, uint32_t &buff_len, bool allow_swap);
			 void PrintEVIOBlockHeader(uint32_t *buff);
		    void PrintStats(void);
			 void DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords=0, const uint32_t *imark=NULL);


		string       source_name;
		bool         is_connected;
		stringstream err_mess;    // last error message
		uint32_t     err_code;    // last error code
		int          ET_STATION_NEVENTS;
		bool         ET_STATION_CREATE_BLOCKING;

		uint64_t     Net_events;
		uint64_t     Nevio_blocks;
		uint64_t     Nevio_events;
		uint64_t     Net_timeouts;
		

#ifdef HAVE_ET
		et_sys_id sys_id;
		et_att_id att_id;
		et_stat_id sta_id;
#endif

		 int VERBOSE;
		bool swap_needed;         // true if block header indicates swapping is needed

		list<pair<uint32_t*, uint32_t> > et_buffs;
		list<pair<uint32_t*, uint32_t> > et_buff_pool;
		
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
};

#endif // _HDET_

