// $Id$
//
//    File: DTOFHit.h
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTOFHit_
#define _DTOFHit_

#include "DObject.h"
#include "DFactory.h"

class DTOFHit:public DObject{
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

