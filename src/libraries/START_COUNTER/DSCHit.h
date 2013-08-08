// $Id$
//
//    File: DSCHit.h
// Created: Wed Feb  7 10:46:20 EST 2007
// Creator: davidl (on Linux megrez.jlab.org 2.6.9-42.0.2.ELsmp x86_64)
//

#ifndef _DSCHit_
#define _DSCHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DSCHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DSCHit);
		
		int sector;	    // sector number 1-24
		float dE;       // Energy loss in GeV
		float t;        // TOF to start counter in ns (best value available)
		float sigma_t;  // uncertainty on t in ns
		bool has_fADC;  // true if this has an fADC hit
		bool has_TDC;   // true if this has an TDC hit

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "sector", "%d", sector);
			AddString(items, "dE", "%3.3f", dE);
			AddString(items, "t", "%3.3f", t);
			AddString(items, "sigma_t", "%d", sigma_t);
			AddString(items, "has_fADC", "%d", has_fADC);
			AddString(items, "has_TDC", "%d", has_TDC);
		}
};

#endif // _DSCHit_

