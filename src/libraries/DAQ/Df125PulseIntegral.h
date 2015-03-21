// $Id$
// $HeadURL$
//
//    File: Df125PulseIntegral.h
// Created: Mon Jul  8 09:17:23 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 x86_64)
//

#ifndef _Df125PulseIntegral_
#define _Df125PulseIntegral_

#include <DAQ/DDAQAddress.h>

class Df125PulseIntegral:public DDAQAddress{
	
	/// Holds pulse integral data for one identified
	/// pulse in one event in one channel of a single
	/// f125 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125PulseIntegral);

		Df125PulseIntegral(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t pulse_number=0, uint32_t quality_factor=0, 
		    uint32_t integral=0, uint32_t pedestal=0, uint32_t nsamples_integral=1, uint32_t nsamples_pedestal=1,bool emulated=false):
		DDAQAddress(rocid, slot, channel, itrigger),pulse_number(pulse_number),quality_factor(quality_factor),
		  integral(integral),pedestal(pedestal),nsamples_integral(nsamples_integral),nsamples_pedestal(nsamples_pedestal),emulated(emulated){}

		uint32_t pulse_number;         ///< from Pulse Integral Data word
		uint32_t quality_factor;       ///< from Pulse Integral Data word
		uint32_t integral;             ///< from Pulse Integral Data word
		uint32_t pedestal;             ///< from Pulse Integral Data word (future)
		uint32_t nsamples_integral;    ///< number of samples used in integral 
		uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
		bool     emulated;             ///< true if made from Window Raw Data
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "quality_factor", "%d", quality_factor);
			AddString(items, "integral", "%d", integral);
			AddString(items, "pedestal", "%d", pedestal);
			AddString(items, "nsamples_integral", "%d", nsamples_integral);
			AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
			AddString(items, "emulated", "%d", emulated);
		}
};

#endif // _Df125PulseIntegral_

