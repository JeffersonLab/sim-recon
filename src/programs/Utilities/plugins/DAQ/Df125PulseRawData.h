// $Id$
// $HeadURL$
//
//    File: Df125PulseRawData.h
// Created: Thu Jun 19 21:08:04 EDT 2014
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df125PulseRawData_
#define _Df125PulseRawData_

#include <DAQ/DDAQAddress.h>

class Df125PulseRawData:public DDAQAddress{
	
	/// Holds pulse raw data for one identified
	/// pulse in one event in one channel of a single
	/// f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df125PulseRawData);
		
		Df125PulseRawData(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t pulse_number=0, uint32_t first_sample_number=0):DDAQAddress(rocid, slot, channel, itrigger),pulse_number(pulse_number),first_sample_number(first_sample_number),invalid_samples(false),overflow(false){}
		
		uint32_t pulse_number;         ///< from Pulse Raw Data Data 1st word
		uint32_t first_sample_number;  ///< from Pulse Raw Data Data 1st word
		vector<uint16_t> samples;      ///< from Pulse Raw Data Data words 2-N (each word contains 2 samples)
		bool invalid_samples;          ///< true if any sample's "not valid" bit set
		bool overflow;                 ///< true if any sample's "overflow" bit set
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "pulse_number", "%d", pulse_number);
			AddString(items, "first_sample_number", "%d", first_sample_number);
			AddString(items, "Nsamples", "%d", samples.size());
			AddString(items, "invalid_samples", "%d", invalid_samples);
			AddString(items, "overflow", "%d", overflow);
		}
};

#endif // _Df125PulseRawData_

