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
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		float MIN_SHOWER_DIST; // not used tight now
                DPi0* makePi0(const TLorentzVector gamma1, const TLorentzVector gamma2, const unsigned int t1, const unsigned int t2); 
};

#endif // _DPi0_factory_

