// $Id$
//
//    File: DBCALHit.h
// Created: Sun Apr  3 10:49:22 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DBCALHit_
#define _DBCALHit_

#include "DFactory.h"

class DBCALHit{
	public:
		HDCLASSDEF(DBCALHit);
		
		enum END_t{
			UPSTREAM,
			DOWNSTREAM
		};
	
		float phim;
		END_t end; /// use BCALHit::UPSTREAM or BCALHit::DOWNSTREAM
		float E;
		float t;
};

#endif // _DBCALHit_

