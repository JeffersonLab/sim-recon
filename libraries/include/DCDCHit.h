// $Id$
//
//    File: DCDCHit.h
// Created: Thu Jun  9 10:22:37 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DCDCHit_
#define _DCDCHit_

#include "DFactory.h"

class DCDCHit:public DObject{
	public:
		HDCLASSDEF(DCDCHit);
		
		float radius;
		float phim;
		float dE;
		float t;
};

#endif // _DCDCHit_

