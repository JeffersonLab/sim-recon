// $Id$
// $HeadURL$
//
//    File: DF1TDCTriggerTime.h
// Created: Sat Nov 24 11:16:07 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DF1TDCTriggerTime_
#define _DF1TDCTriggerTime_

#include <DAQ/DDAQAddress.h>

class DF1TDCTriggerTime:public DDAQAddress{
	
	/// Holds trigger time data for one event in
	/// a single F1 TDC module.
	
	public:
		JOBJECT_PUBLIC(DF1TDCTriggerTime);

		DF1TDCTriggerTime():DDAQAddress(0, 0, 0, 0),time(0){}
		DF1TDCTriggerTime(uint32_t rocid, uint32_t slot, uint32_t itrigger, uint64_t time):DDAQAddress(rocid, slot, 0, itrigger),time(time){}
		
		uint64_t time;           // from Trigger Time words
		
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "time", "%ld", time);
		}

};

#endif // _DF1TDCTriggerTime_

