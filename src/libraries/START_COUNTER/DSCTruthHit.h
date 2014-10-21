// $Id$
//
//    File: DSCTruthHit.h
// Created: Wed Feb  7 10:53:46 EST 2007
// Creator: davidl (on Linux megrez.jlab.org 2.6.9-42.0.2.ELsmp x86_64)
//

#ifndef _DSCTruthHit_
#define _DSCTruthHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DSCTruthHit:public JObject{
	public:
		JOBJECT_PUBLIC(DSCTruthHit);
		
		float dEdx;
		bool primary;
		int track;
		int itrack;
                int ptype;
		float r;
		float phi;
		float z;
		float t;
		int sector;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "track", "%d", track);
			AddString(items, "itrack", "%d", itrack);
			AddString(items, "primary", "%d", primary);
                        AddString(items, "ptype", "%d", ptype);
			AddString(items, "dEdx(MeV/cm)", "%1.3f", dEdx*1.0E3);
			AddString(items, "t", "%3.2f", t);
			AddString(items, "r", "%3.1f", r);
			AddString(items, "phi", "%1.3f", phi);
			AddString(items, "z", "%3.1f", z);
                        AddString(items, "sector", "%d", sector);
		}
};

#endif // _DSCTruthHit_

