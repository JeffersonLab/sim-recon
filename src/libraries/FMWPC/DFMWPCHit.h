// $Id$
//
//    File: DFMWPCHit.h
// Created: Tue Jun 16 07:07:54 EDT 2015
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DFMWPCHit_
#define _DFMWPCHit_

#include <JANA/jerror.h>

class DFMWPCHit:public JObject{
	public:
		JOBJECT_PUBLIC(DFMWPCHit);

		int layer;   // 1-8
		int wire;    // 1-144
		float dE;    // energy in GeV
		float t;     // time in ns

		void toStrings(vector<pair<string, string> >&items) const {
			AddString(items, "layer", "%d", layer);
			AddString(items, "wire", "%d", wire);
			AddString(items, "dE(keV)", "%3.1f", dE*1.0E6);
			AddString(items, "t", "%3.3f", t);
		}

};

#endif // _DFMWPCHit_

