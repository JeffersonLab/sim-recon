// $Id$
//
//    File: DF1TDCTriggerTime.h
// Created: Sat Nov 24 11:16:07 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DF1TDCTriggerTime_
#define _DF1TDCTriggerTime_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>

class DF1TDCTriggerTime:public jana::JObject{
	
	/// Holds trigger time data for one event in
	/// a single F1 TDC module.
	
	public:
		JOBJECT_PUBLIC(DF1TDCTriggerTime);

		DF1TDCTriggerTime():rocid(0),slot(0),itrigger(0),time(0){}
		DF1TDCTriggerTime(uint32_t rocid, uint32_t slot, uint32_t itrigger, uint64_t time):rocid(rocid),slot(slot),itrigger(itrigger),time(time){}
		
		uint32_t rocid;          // from EVIO header (crate number)
		uint32_t slot;           // from Block Header
		uint32_t itrigger;       // from Event Header
		uint64_t time;           // from Trigger Time words
		
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid", "%d", rocid);
			AddString(items, "slot", "%d", slot);
			AddString(items, "itrigger", "%d", itrigger);
			AddString(items, "time", "%ld", time);
		}

};

#endif // _DF1TDCTriggerTime_

