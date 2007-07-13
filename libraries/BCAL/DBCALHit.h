// $Id$
//
//    File: DBCALHit.h
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DBCALHit_
#define _DBCALHit_

#include "BCAL/DBCALGeometry.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DBCALHit:public JObject{
	public:
		HDCLASSDEF(DBCALHit);
		
		int module;
		int layer;
		int sector;
		DBCALGeometry::End end;
		float E;
		float t;
};

#endif // _DBCALHit_

