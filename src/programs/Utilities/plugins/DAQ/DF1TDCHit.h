// $Id$
//
//    File: DF1TDCHit.h
// Created: Fri Aug 10 12:02:49 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DF1TDCHit_
#define _DF1TDCHit_

#include <DAQ/DDAQAddress.h>

class DF1TDCHit:public DDAQAddress{
	
	/// Holds single hit from a F1TDC module
	
	public:
		JOBJECT_PUBLIC(DF1TDCHit);
		
		DF1TDCHit(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t trig_time=0, uint32_t time=0):DDAQAddress(rocid, slot, channel, itrigger),trig_time(trig_time),time(time){}
		
		uint32_t trig_time;            // from data word
		uint32_t time;                 // from data word
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "trig_time", "%d", trig_time);
			AddString(items, "time", "%d", time);
		}
};

#endif // _DF1TDCHit_

