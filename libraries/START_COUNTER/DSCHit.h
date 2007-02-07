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

class DSCHit:public JObject{
	public:
		HDCLASSDEF(DSCHit);
		
		float dE;		// Energy loss in GeV
		float t;			// TOF to start counter in ns
		int sector;		// sector number 1-24
};

#endif // _DSCHit_

