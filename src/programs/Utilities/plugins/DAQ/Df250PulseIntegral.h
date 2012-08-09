// $Id$
//
//    File: Df250PulseIntegral.h
// Created: Tue Aug  7 15:24:50 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df250PulseIntegral_
#define _Df250PulseIntegral_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>

class Df250PulseIntegral:public jana::JObject{
	
	/// Holds pulse integral data for one identified
	/// pulse in one event in one channel of a single
	/// f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250PulseIntegral);

		Df250PulseIntegral():rocid(0),slot(0),channel(0),pulse_number(0),quality_factor(0),integral(0){}
		
		uint32_t rocid;                // from EVIO header (crate number)
		uint32_t slot;                 // from Block Header
		uint32_t channel;              // from Pulse Integral Data word
		uint32_t pulse_number;         // from Pulse Integral Data word
		uint32_t quality_factor;       // from Pulse Integral Data word
		uint32_t integral;             // from Pulse Integral Data word
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid", "%d", rocid);
			AddString(items, "slot", "%d", slot);
			AddString(items, "channel", "%d", channel);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "quality_factor", "%d", quality_factor);
			AddString(items, "integral", "%d", integral);
		}
};

#endif // _Df250PulseIntegral_

