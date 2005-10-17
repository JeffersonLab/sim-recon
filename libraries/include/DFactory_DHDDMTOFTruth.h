// $Id$
//
//    File: DFactory_DHDDMTOFTruth.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFactory_DHDDMTOFTruth_
#define _DFactory_DHDDMTOFTruth_

#include "DFactory.h"
#include "DHDDMTOFTruth.h"

class DFactory_DHDDMTOFTruth:public DFactory<DHDDMTOFTruth>{
    public:
        DFactory_DHDDMTOFTruth(){};
        ~DFactory_DHDDMTOFTruth(){};

        derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		
        const string toString(void);


    private:
        derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DHDDMTOFTruth_

