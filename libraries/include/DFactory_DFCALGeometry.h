// $Id$
//
//    File: DFactory_DFCALGeometry.h
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFactory_DFCALGeometry_
#define _DFactory_DFCALGeometry_

#include "DFactory.h"
#include "DFCALGeometry.h"

class DFactory_DFCALGeometry : public DFactory<DFCALGeometry> {

public:
	
	DFactory_DFCALGeometry() { flags = PERSISTANT; }
	~DFactory_DFCALGeometry(){}
	const string toString(void);

private:
	
	derror_t evnt(DEventLoop *loop, int eventnumber);	
};

#endif // _DFactory_DFCALGeometry_

