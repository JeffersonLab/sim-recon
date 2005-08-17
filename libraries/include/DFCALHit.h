// $Id$
//
//    File: DFCALHit.h
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFCALHit_
#define _DFCALHit_

#include "DObject.h"
#include "DFactory.h"

class DFCALHit:public DObject{
	public:
		HDCLASSDEF(DFCALHit);
		
		float x;
		float y;
		float E;
		float t;
};

#endif // _DFCALHit_

