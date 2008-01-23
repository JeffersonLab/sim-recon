// $Id$
//
//    File: DFCALMCResponse_factory.h
// Created: Wed Aug 24 10:05:31 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#ifndef _DFCALMCResponse_factory_
#define _DFCALMCResponse_factory_

#include "JANA/JFactory.h"
#include "DFCALMCResponse.h"

#include "DRandom.h"

class DFCALMCResponse_factory:public JFactory<DFCALMCResponse>{

public:
	
	DFCALMCResponse_factory();
	~DFCALMCResponse_factory(){};

	const string toString(void);

private:
		
		jerror_t evnt(JEventLoop *loop, int eventnumber);	
 
    float m_photStatCoef; // photon statistis contribution to energy resolution
    float m_blockThreshold;

    DRandom m_randomGen;

};

#endif // _DFCALMCResponse_factory_

