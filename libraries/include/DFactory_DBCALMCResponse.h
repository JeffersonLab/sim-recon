// $Id$
//
//    File: DFactory_DBCALMCResponse.h
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DFactory_DBCALMCResponse_
#define _DFactory_DBCALMCResponse_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DBCALMCResponse.h"

class DFactory_DBCALMCResponse:public DFactory<DBCALMCResponse>{
	public:
		DFactory_DBCALMCResponse(){};
		~DFactory_DBCALMCResponse(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DBCALMCResponse_

