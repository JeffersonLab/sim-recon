// $Id$
//
//    File: DHDDMTOFTruth_factory.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DHDDMTOFTruth_factory_
#define _DHDDMTOFTruth_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"
#include "DHDDMTOFTruth.h"

class DHDDMTOFTruth_factory:public JFactory<DHDDMTOFTruth>{
    public:
        DHDDMTOFTruth_factory(){};
        ~DHDDMTOFTruth_factory(){};

        jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		
        const string toString(void);


    private:
        jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DHDDMTOFTruth_factory_

