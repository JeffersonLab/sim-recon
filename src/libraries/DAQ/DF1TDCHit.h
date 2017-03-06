// $Id$
// $HeadURL$
//
//    File: DF1TDCHit.h
// Created: Fri Aug 10 12:02:49 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DF1TDCHit_
#define _DF1TDCHit_

#include <DAQ/DDAQAddress.h>
#include "DAQ/DModuleType.h"

class DF1TDCHit:public DDAQAddress{
	
	/// Holds single hit from a F1TDC module
	
	public:
		JOBJECT_PUBLIC(DF1TDCHit);
		
		DF1TDCHit(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t trig_time=0, uint32_t time=0, uint32_t data_word=0, MODULE_TYPE locModType = DModuleType::UNKNOWN) :
		DDAQAddress(rocid, slot, channel, itrigger), trig_time(trig_time), time(time), data_word(data_word), modtype(locModType) {
			res_status                  = (data_word>>26) & 0x1;
			output_fifo_overflow_status = (data_word>>25) & 0x1; 
			hit_fifo_overflow_status    = (data_word>>24) & 0x1; 
		}

		uint32_t trig_time;                    // from data word
		uint32_t time;                         // from data word
		bool     res_status;                   // 1 = locked    0 = unlocked
		bool     output_fifo_overflow_status;  // 1 = overflow  0 = normal
		bool     hit_fifo_overflow_status;     // 1 = overflow  0 = normal
		uint32_t data_word;                    // full data word (bits 24-26 contain some status info)
		MODULE_TYPE modtype;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "trig_time", "%d", trig_time);
			AddString(items, "time", "%d", time);
			AddString(items, "data_word", "0x%08x", data_word);
		}
};

#endif // _DF1TDCHit_

