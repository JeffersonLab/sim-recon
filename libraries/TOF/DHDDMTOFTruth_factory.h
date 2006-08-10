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

		
        const string toString(void);


    private:
};

#endif // _DHDDMTOFTruth_factory_

