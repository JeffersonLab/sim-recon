// $Id$
//
//    File: DTOFHit.h
// Created: Sun Apr  3 10:31:26 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTOFHit_
#define _DTOFHit_

#include "DFactory.h"

class DTOFHit{
	public:
		HDCLASSDEF(DTOFHit);
		
		float x;
		float y;
		float dE;
		float t;
		int orientation;	///< 0=vertical  1=horizontal
		int end;				///< 0=left/top 1=right/bottom
};

#endif // _DTOFHit_

