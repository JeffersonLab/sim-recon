// $Id$
//
//    File: DFactory_DTOFGeometry.h
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DFactory_DTOFGeometry_
#define _DFactory_DTOFGeometry_

#include "DFactory.h"
#include "DTOFGeometry.h"

class DFactory_DTOFGeometry:public DFactory<DTOFGeometry>{
	public:
		DFactory_DTOFGeometry(){};
		~DFactory_DTOFGeometry(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTOFGeometry_

