// $Id$
//
//    File: DFCALHit.h
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFCALHit_
#define _DFCALHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DFCALHit:public JObject{
	
	public:
		HDCLASSDEF(DFCALHit);
		DFCALHit(){}
		
	DFCALHit( oid_t id, 
			  float x, float y,
			  float E, float t ) :
		JObject( id ),
		x( x ),
		y( y ),
		E( E ),
		t( t ) {}
		
		int column;
		int row;
		float x;
		float y;
		float E;
		float t;
};

#endif // _DFCALHit_

