// $Id$
//
//    File: DUPVHit.h
// Created: Thu Jun  9 10:01:38 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DUPVHit_
#define _DUPVHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DUPVHit:public JObject{
	public:
		HDCLASSDEF(DUPVHit);
		
		enum UPV_side_t{
			UPV_LEFT = 0,
			UPV_RIGHT = 1
		};
		
		int layer;
		int row;
		float E;		// GeV
		float t;		// ns
		int side;	// 0=left 1=right
		
};

#endif // _DUPVHit_

