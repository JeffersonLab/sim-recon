// $Id$
//
//    File: DTOFTruth_factory.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFTruth_factory_
#define _DTOFTruth_factory_

#include "JANA/JFactory.h"

#include "DTOFTruth.h"

class DTOFTruth_factory:public JFactory<DTOFTruth>{
    public:
        DTOFTruth_factory(){};
        ~DTOFTruth_factory(){};

		
        const string toString(void);


    private:
};

#endif // _DTOFTruth_factory_

