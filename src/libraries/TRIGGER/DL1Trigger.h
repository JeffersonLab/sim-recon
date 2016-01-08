// $Id$
//
//    File: DL1Trigger.h
// Created: Fri Jan  8 10:57:58 EST 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DL1Trigger_
#define _DL1Trigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DL1Trigger:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DL1Trigger);
		
		DL1Trigger():gtp_latch(0),fp_latch(0){}
		
		uint32_t gtp_latch; // latch word of active triggers from Global Trigger Processor
		uint32_t fp_latch;  // latch word of active triggers from Front Panel
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "gtp_latch", "%04x", gtp_latch);
			AddString(items, "fp_latch",  "%04x", fp_latch);
		}
		
};

#endif // _DL1Trigger_

