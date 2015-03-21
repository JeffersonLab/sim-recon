// $Id$
// $HeadURL$
//
//    File: Df125TriggerTime.h
// Created: Mon Jul  8 11:58:08 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df125TriggerTime_
#define _Df125TriggerTime_

#include <DAQ/DDAQAddress.h>

class Df125TriggerTime:public DDAQAddress{
	
	/// Holds trigger time data for one event in
	/// a single f125 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125TriggerTime);

		Df125TriggerTime():DDAQAddress(0, 0, 0),itrigger(0),time(0){}
		Df125TriggerTime(uint32_t rocid, uint32_t slot, uint32_t itrigger, uint64_t time):DDAQAddress(rocid, slot, 0, itrigger),time(time){}
		
		uint32_t itrigger;       // from Event Header
		uint64_t time;           // from Trigger Time words
		
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "time", "%ld", time);
		}

};

#endif // _Df125TriggerTime_

