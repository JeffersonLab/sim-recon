// $Id$
//
//    File: DFCALGeometry_factory.h
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALGeometry_factory_
#define _DFCALGeometry_factory_

#include "JANA/JFactory.h"
#include "DFCALGeometry.h"

class DFCALGeometry_factory : public JFactory<DFCALGeometry> {

public:
	
	DFCALGeometry_factory() {}
	~DFCALGeometry_factory(){}

private:
	
	jerror_t brun(JEventLoop *loop, int32_t runnumber);	
	jerror_t erun(void);	
};

#endif // _DFCALGeometry_factory_

