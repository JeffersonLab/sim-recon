// $Id$
//
//    File: DFactory_DFCALMCResponse.h
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFactory_DFCALMCResponse_
#define _DFactory_DFCALMCResponse_

#include "DFactory.h"
#include "DFCALMCResponse.h"

class DFactory_DFCALMCResponse:public DFactory<DFCALMCResponse>{

public:
	
	DFactory_DFCALMCResponse(){};
	~DFactory_DFCALMCResponse(){};

	const string toString(void);

private:
		
		derror_t evnt(DEventLoop *loop, int eventnumber);	
};

#endif // _DFactory_DFCALMCResponse_

