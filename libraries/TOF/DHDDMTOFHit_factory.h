// $Id$
//
//    File: DHDDMTOFHit_factory.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DHDDMTOFHit_factory_
#define _DHDDMTOFHit_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"

#include "DHDDMTOFHit.h"

class DHDDMTOFHit_factory:public JFactory<DHDDMTOFHit>{
    public:

        DHDDMTOFHit_factory(){};
        ~DHDDMTOFHit_factory(){};
		
        const string toString(void);

    private:
        jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DHDDMTOFHit_factory_

