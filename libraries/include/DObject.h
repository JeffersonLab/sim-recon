// $Id$
//
//    File: DObject.h
// Created: Wed Aug 17 10:57:09 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DObject_
#define _DObject_

#include "derror.h"


typedef int identifier_t;


class DObject{

	/// The DObject class is a base class for all data classes.
	/// In other words, all classes which are produced by DANA
	/// factories.

	public:
	
		DObject();
   	    DObject( identifier_t aId ) : id( aId ) {}

		virtual ~DObject(){}

		identifier_t id;
		
};

#endif // _DObject_

