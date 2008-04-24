//
//    File: DPi0_factory.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DPi0_factory_
#define _DPi0_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DPi0.h"

class DPi0_factory:public JFactory<DPi0>{
	public:
		DPi0_factory();
		~DPi0_factory(){};
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		float MIN_SHOWER_DIST; // not used right now
                //DPi0* makePi0(const DLorentzVector& gamma1, const DLorentzVector& gamma2, unsigned int t1, unsigned int t2, oid_t id1, oid_t id2); 
                DPi0* makePi0(const DPhoton* photon1, const DPhoton* photon2); 
};

		// void makePi0fromID( oid_t child1, oid_t child2);

#endif // _DPi0_factory_

