// $Id$
//
//    File: Df250WindowSum.h
// Created: Tue Aug  7 15:24:31 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _Df250WindowSum_
#define _Df250WindowSum_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>

class Df250WindowSum:public jana::JObject{
	
	/// Holds window sum data for one event in
	/// one channel of a single f250 Flash ADC module.
	
	public:
		JOBJECT_PUBLIC(Df250WindowSum);
	
		Df250WindowSum():rocid(0),slot(0),channel(0),sum(0),overflow(false){}

		uint32_t rocid;          // from EVIO header (crate number)
		uint32_t slot;           // from Block Header
		uint32_t channel;        // from Window Sum Data word
		uint32_t sum;            // from Window Sum Data word
		bool overflow;           // true if "overflow" bit set
	
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "rocid", "%d", rocid);
			AddString(items, "slot", "%d", slot);
			AddString(items, "channel", "%d", channel);
			AddString(items, "sum", "%d", sum);
			AddString(items, "overflow", "%d", overflow);
		}

};

#endif // _Df250WindowSum_

