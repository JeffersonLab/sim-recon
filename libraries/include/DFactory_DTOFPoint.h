// $Id$
//
//    File: DFactory_DTOFPoint.h
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFactory_DTOFPoint_
#define _DFactory_DTOFPoint_

#include "DFactory.h"
#include "DTOFPoint.h"

class DFactory_DTOFPoint:public DFactory<DTOFPoint>{
	public:
		DFactory_DTOFPoint(){};
		~DFactory_DTOFPoint(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTOFPoint_

