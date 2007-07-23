/*
 *  DMCBCALHit_factory.h
 *  Hall D
 *
 *  Created by Matthew Shepherd on 7/23/07.
 *
 */


#ifndef _DMCBCALHit_factory_
#define _DMCBCALHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DHDDMBCALHit.h"

class DMCBCALHit_factory : public JFactory< DHDDMBCALHit > { 
    
public:
    
    DMCBCALHit_factory();
    ~DMCBCALHit_factory(){};
    
    const string toString( void );
    
};

#endif
