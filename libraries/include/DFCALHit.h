// $Id$
//
//    File: DFCALHit.h
// Created: Sun Apr  3 10:41:55 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFCALHit_
#define _DFCALHit_

#include "DFactory.h"

class DFCALHit{
	public:
		HDCLASSDEF(DFCALHit);
		
		float x;
		float y;
		float E;
		float t;
};

#endif // _DFCALHit_

