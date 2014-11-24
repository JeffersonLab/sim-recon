// $Id$
// $HeadURL$
//
//    File: Df250WindowSum.h
// Created: Tue Aug  7 15:24:31 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df250WindowSum_
#define _Df250WindowSum_

#include <DAQ/DDAQAddress.h>

class Df250WindowSum:public DDAQAddress{
	
	/// Holds window sum data for one event in
	/// one channel of a single f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250WindowSum);
	
		Df250WindowSum(uint32_t rocid=0, uint32_t slot=0, uint32_t channel=0, uint32_t itrigger=0, uint32_t sum=0, bool overflow=false):DDAQAddress(rocid, slot, channel, itrigger),sum(sum),overflow(overflow){}

		uint32_t sum;            // from Window Sum Data word
		bool overflow;           // true if "overflow" bit set
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			DDAQAddress::toStrings(items);
			AddString(items, "sum", "%d", sum);
			AddString(items, "overflow", "%d", overflow);
		}

};

#endif // _Df250WindowSum_

