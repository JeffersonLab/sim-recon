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

#include "DCDCWire.h"

class DCDCTrackHit:public JObject{
	public:
		HDCLASSDEF(DCDCTrackHit);

		const DCDCWire *wire;	// stereo angle of tube in radians
		float tdrift;				// Drift time of hit in ns
		float dist;					// Measured DOCA in cm
};

#endif // _DCDCTrackHit_

