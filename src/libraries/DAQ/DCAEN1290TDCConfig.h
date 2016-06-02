// $Id$
// $HeadURL$
//
//    File: DCAEN1290TDCConfig.h
// Created: Sun Sep  7 15:52:50 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _DCAEN1290TDCConfig_
#define _DCAEN1290TDCConfig_

#include <DAQ/DDAQConfig.h>

class DCAEN1290TDCConfig:public DDAQConfig{
	public:
		JOBJECT_PUBLIC(DCAEN1290TDCConfig);
		
		DCAEN1290TDCConfig(uint32_t rocid, uint32_t slot_mask):DDAQConfig(rocid,slot_mask),WINWIDTH(0xFFFF),WINOFFSET(0xFFFF){}
		DCAEN1290TDCConfig(const DCAEN1290TDCConfig *c):DDAQConfig(c->rocid,c->slot_mask),WINWIDTH(c->WINWIDTH),WINOFFSET(c->WINOFFSET){}
		
		uint16_t WINWIDTH;
		uint16_t WINOFFSET;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQConfig::toStrings(items);
			AddString(items, "WINWIDTH"  , "%d", WINWIDTH);
			AddString(items, "WINOFFSET" , "%d", WINOFFSET);
		}
		
};

#endif // _DCAEN1290TDCConfig_

