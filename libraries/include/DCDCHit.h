// $Id$
//
//    File: DCDCHit.h
// Created: Sun Apr  3 10:46:28 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DCDCHit_
#define _DCDCHit_

#include "DFactory.h"

class DCDCHit{
	public:
		HDCLASSDEF(DCDCHit);
		
		float radius;
		float phim;
		float dE;
		float t;
};

#endif // _DCDCHit_

