// $Id$
//
//    File: DFactory_DHDDMTOFHit.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFactory_DHDDMTOFHit_
#define _DFactory_DHDDMTOFHit_

#include "DFactory.h"
#include "DHDDMTOFHit.h"

class DFactory_DHDDMTOFHit:public DFactory<DHDDMTOFHit>{
    public:
        DFactory_DHDDMTOFHit(){};
        ~DFactory_DHDDMTOFHit(){};

        derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		
        const string toString(void);


    private:
        derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DHDDMTOFHit_

