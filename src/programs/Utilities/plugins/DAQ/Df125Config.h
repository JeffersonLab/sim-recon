// $Id$
// $HeadURL$
//
//    File: Df125Config.h
// Created: Sun Sep  7 15:52:50 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _Df125Config_
#define _Df125Config_

#include <DAQ/DDAQConfig.h>

class Df125Config:public DDAQConfig{
	public:
		JOBJECT_PUBLIC(Df125Config);

		Df125Config(uint32_t rocid, uint32_t slot_mask):DDAQConfig(rocid,slot_mask),NSA(0xFFFF),NSB(0xFFFF),NSA_NSB(0xFFFF),NPED(0xFFFF),WINWIDTH(0xFFFF){}
		
		uint16_t NSA;      // Num. samples before threshold crossing sample
		uint16_t NSB;      // Num. samples after  threshold crossing sample
		uint16_t NSA_NSB;  // NSA+NSB = total number of samples in integration window
		uint16_t NPED;     // Number of samples used to determine pedestal
		uint16_t WINWIDTH; // maximum integration window size (in samples)
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQConfig::toStrings(items);
			AddString(items, "NSA"      , "%d", NSA);
			AddString(items, "NSB"      , "%d", NSB);
			AddString(items, "NSA_NSB"  , "%d", NSA_NSB);
			AddString(items, "NPED"     , "%d", NPED);
			AddString(items, "WINWIDTH" , "%d", WINWIDTH);
		}
		
};

#endif // _Df125Config_

