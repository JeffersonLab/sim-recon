// $Id$
//
//    File: Df250PulseTime.h
// Created: Tue Aug  7 15:25:03 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df250PulseTime_
#define _Df250PulseTime_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>

class Df250PulseTime:public jana::JObject{
	
	/// Holds pulse time for one identified
	/// pulse in one event in one channel of a single
	/// f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250PulseTime);
		
		Df250PulseTime(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t pulse_number=0, uint32_t quality_factor=0, uint32_t time=0):rocid(rocid),slot(slot),channel(channel),pulse_number(pulse_number),quality_factor(quality_factor),time(time){}
		
		uint32_t rocid;                // from EVIO header (crate number)
		uint32_t slot;                 // from Block Header
		uint32_t channel;              // from Pulse Time Data word
		uint32_t pulse_number;         // from Pulse Time Data word
		uint32_t quality_factor;       // from Pulse Time Data word
		uint32_t time;                 // from Pulse Time Data word
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid", "%d", rocid);
			AddString(items, "slot", "%d", slot);
			AddString(items, "channel", "%d", channel);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "quality_factor", "%d", quality_factor);
			AddString(items, "time", "%d", time);
		}
};

#endif // _Df250PulseTime_

