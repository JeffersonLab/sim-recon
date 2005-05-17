// $Id$
//
//    File: DFactory_DFCALShowers.h
// Created: Tue May 17 09:47:34 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFactory_DFCALShowers_
#define _DFactory_DFCALShowers_

#include "DFactory.h"
#include "DFCALShowers.h"

class DFactory_DFCALShowers:public DFactory<DFCALShowers>{
	public:
		DFactory_DFCALShowers(){};
		~DFactory_DFCALShowers(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DFCALShowers_

