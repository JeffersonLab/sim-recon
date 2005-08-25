// $Id$
//
//    File: DTrackHit.h
// Created: Tue Aug 23 05:00:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackHit_
#define _DTrackHit_

#include "DObject.h"
#include "DFactory.h"
#include "GlueX.h"

class DTrackHit:public DObject{
	public:
		HDCLASSDEF(DTrackHit);

		float x,y,z,r,phi;
		DetectorSystem_t system;
};

#endif // _DTrackHit_

