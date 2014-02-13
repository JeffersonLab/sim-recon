// $Id$
//
//    File: DBCALGeometry_factory.h
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALGeometry_factory_
#define _DBCALGeometry_factory_

#include <JANA/JFactory.h>
using namespace jana;

#include "DBCALGeometry.h"

class DBCALGeometry_factory:public JFactory<DBCALGeometry>{
	public:
		DBCALGeometry_factory(){};
		~DBCALGeometry_factory(){};


	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DBCALGeometry_factory_

