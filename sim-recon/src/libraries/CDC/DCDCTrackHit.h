// $Id$
//
//    File: DCDCTrackHit.h
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DCDCTrackHit_
#define _DCDCTrackHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#include "DCDCWire.h"

class DCDCTrackHit:public JObject{
	public:
		JOBJECT_PUBLIC(DCDCTrackHit);

		const DCDCWire *wire;	// DCDCWire structure for this hit
		bool is_stereo; // true if this is stereo wire
		float tdrift;				// Drift time of hit in ns
		float dist;					// Measured DOCA in cm
		float dE; // Energy deposition in GeV

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "ring", "%d", wire->ring);
			AddString(items, "straw", "%d", wire->straw);
			AddString(items, "x(cm)", "%3.1f", wire->origin.x());
			AddString(items, "y(cm)", "%3.1f", wire->origin.y());
			AddString(items, "stereo(rad)", "%1.4f", wire->stereo);
			AddString(items, "tdrift(ns)", "%3.1f", tdrift);
			AddString(items, "dist(cm)", "%1.3f", dist);
		}
};

#endif // _DCDCTrackHit_

