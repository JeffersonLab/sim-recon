// $Id$
//
//    File: DObject.h
// Created: Wed Aug 17 10:57:09 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DObject_
#define _DObject_

#include "derror.h"


typedef unsigned long oid_t;


class DObject{

	/// The DObject class is a base class for all data classes.
	/// In other words, all classes which are produced by DANA
	/// factories.

	public:
	
		DObject(){id = (oid_t)this;}
		DObject( oid_t aId ) : id( aId ) {}

		virtual ~DObject(){}
		
		template<typename T>
		bool IsA(const T *t){return dynamic_cast<const T*>(this)!=0L;}

		oid_t id;
		
};

#endif // _DObject_

