// $Id$
//
//    File: DFactory_DHDDMForwardShower.h
// Created: Mon Aug 29 15:14:08 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFactory_DHDDMForwardShower_
#define _DFactory_DHDDMForwardShower_

#include "DFactory.h"
#include "DHDDMForwardShower.h"

class DFactory_DHDDMForwardShower : public DFactory<DHDDMForwardShower>{

public:
	
	DFactory_DHDDMForwardShower(){};
	~DFactory_DHDDMForwardShower(){};

	derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		
	const string toString(void);

private:
		
	derror_t evnt(DEventLoop *loop, int eventnumber);
};

#endif // _DFactory_DHDDMForwardShower_

