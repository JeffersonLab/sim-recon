// $Id$
//
//    File: Df125PulseTime.h
// Created: Mon Jul  8 12:01:57 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df125PulseTime_
#define _Df125PulseTime_

#include <DAQ/DDAQAddress.h>

class Df125PulseTime:public DDAQAddress{
	
	/// Holds pulse time for one identified
	/// pulse in one event in one channel of a single
	/// f125 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125PulseTime);
		
		Df125PulseTime(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t pulse_number=0, uint32_t time=0):DDAQAddress(rocid, slot, channel, itrigger),pulse_number(pulse_number),time(time){}
		
		uint32_t pulse_number;         // from Pulse Time Data word
		uint32_t time;                 // from Pulse Time Data word
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "time", "%d", time);
		}
};

#endif // _Df125PulseTime_

