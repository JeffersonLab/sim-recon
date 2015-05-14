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
		
		int sector;	    // sector number 1-30
		float dE;       // Energy loss in GeV
		float t;        // best time (walk-corrected tdc)
		float t_TDC;   // time from TDC, no walk correction
		float t_fADC; // time from fADC
		float pulse_height; // amplitude of pulse (used in time-walk corrections)
		bool has_fADC; 
		bool has_TDC;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "sector", "%d", sector);
			AddString(items, "dE", "%3.3f", dE);
			AddString(items, "t", "%3.3f", t);
			AddString(items, "t_TDC","%3.3f", t_TDC);
			AddString(items, "t_fADC", "%3.3f", t_fADC);
			AddString(items, "has_fADC", "%d", (int)has_fADC);
			AddString(items, "has_TDC", "%d", (int)has_TDC);
		}
};

#endif // _DSCHit_

