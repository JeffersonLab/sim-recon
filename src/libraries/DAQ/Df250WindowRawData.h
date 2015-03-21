// $Id$
// $HeadURL$
//
//    File: Df250WindowRawData.h
// Created: Tue Aug  7 15:24:27 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df250WindowRawData_
#define _Df250WindowRawData_

#include <DAQ/DDAQAddress.h>

class Df250WindowRawData:public DDAQAddress{

	/// Holds window raw data samples for one event in
	/// one channel of a single f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250WindowRawData);
	
		Df250WindowRawData(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0):DDAQAddress(rocid, slot, channel, itrigger),invalid_samples(false),overflow(false){}
	
		vector<uint16_t> samples;// from Window Raw Data words 2-N (each word contains 2 samples)
		bool invalid_samples;    // true if any sample's "not valid" bit set
		bool overflow;           // true if any sample's "overflow" bit set

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "Nsamples", "%d", samples.size());
			AddString(items, "invalid_samples", "%d", invalid_samples);
			AddString(items, "overflow", "%d", overflow);
		}

};

#endif // _Df250WindowRawData_

