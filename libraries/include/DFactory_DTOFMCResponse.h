// $Id$
//
//    File: DFactory_DTOFMCResponse.h
// Created: Mon Aug 15 11:33:45 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#ifndef _DFactory_DTOFMCResponse_
#define _DFactory_DTOFMCResponse_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DTOFMCResponse.h"

class DFactory_DTOFMCResponse:public DFactory<DTOFMCResponse>{
	public:
		DFactory_DTOFMCResponse(){};
		~DFactory_DTOFMCResponse(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTOFMCResponse_

