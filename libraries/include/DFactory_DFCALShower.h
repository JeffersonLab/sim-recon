// $Id$
//
//    File: DFactory_DFCALShower.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFactory_DFCALShower_
#define _DFactory_DFCALShower_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DFCALShower.h"

class DFactory_DFCALShower:public DFactory<DFCALShower>{
	public:
		DFactory_DFCALShower(){};
		~DFactory_DFCALShower(){};
		const string toString(void);
	
	private:
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DFCALShower_

