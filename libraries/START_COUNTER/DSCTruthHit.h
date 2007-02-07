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
		HDCLASSDEF(DSCTruthHit);
		
		float dEdx;
		bool primary;
		int track;
		float r;
		float phi;
		float z;
		float t;
};

#endif // _DSCTruthHit_

