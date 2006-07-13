// $Id$
//
//    File: DBCALHit.h
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DBCALHit_
#define _DBCALHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DBCALHit:public JObject{
	public:
		HDCLASSDEF(DBCALHit);
		
		enum END_t{
			UPSTREAM,
			DOWNSTREAM
		};
		
		int module;
		int layer;
		int sector;
		END_t end; /// use BCALHit::UPSTREAM or BCALHit::DOWNSTREAM
		float E;
		float t;
};

#endif // _DBCALHit_

