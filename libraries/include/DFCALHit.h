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
		
	DFCALHit( identifier_t id, 
			  float x, float y,
			  float E, float t ) :
		DObject( id ),
		x( x ),
		y( y ),
		E( E ),
		t( t ) {}
		
		float x;
		float y;
		float E;
		float t;
};

#endif // _DFCALHit_

