// $Id$
// $HeadURL$
//
//    File: Df125PulsePedestal.h
// Created: Mon Jul 28 09:44:35 EDT 2014
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df125PulsePedestal_
#define _Df125PulsePedestal_

#include <DAQ/DDAQAddress.h>

class Df125PulsePedestal:public DDAQAddress{
	
	/// Holds pulse time for one identified
	/// pulse in one event in one channel of a single
	/// f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125PulsePedestal);
		
		Df125PulsePedestal(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t pulse_number=0, uint32_t pedestal=0, uint32_t pulse_peak=0, uint32_t nsamples=0, bool emulated=false):DDAQAddress(rocid, slot, channel, itrigger),pulse_number(pulse_number),pedestal(pedestal),pulse_peak(pulse_peak),nsamples(nsamples),emulated(emulated){}
		
		uint32_t pulse_number;   ///< from Pulse Pedestal Data word
		uint32_t pedestal;       ///< from Pulse Pedestal Data word
		uint32_t pulse_peak;     ///< from Pulse Pedestal Data word
		uint32_t nsamples;       ///< number of samples used in pedestal
		bool     emulated;       ///< true if made from Window Raw Data
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "pedestal", "%d", pedestal);
			AddString(items, "pulse_peak", "%d", pulse_peak);
			AddString(items, "nsamples", "%d", nsamples);
			AddString(items, "emulated", "%d", emulated);
		}
};

#endif // _Df125PulsePedestal_

