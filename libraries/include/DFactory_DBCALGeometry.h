// $Id$
//
//    File: DFactory_DBCALGeometry.h
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DFactory_DBCALGeometry_
#define _DFactory_DBCALGeometry_

#include "DFactory.h"
#include "DBCALGeometry.h"

class DFactory_DBCALGeometry:public DFactory<DBCALGeometry>{
	public:
		DFactory_DBCALGeometry(){};
		~DFactory_DBCALGeometry(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DBCALGeometry_

