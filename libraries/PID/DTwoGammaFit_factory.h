//
//    File: DTwoGammaFit_factory.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_factory_
#define _DTwoGammaFit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTwoGammaFit.h"

class DTwoGammaFit_factory:public JFactory<DTwoGammaFit>{
	public:
		DTwoGammaFit_factory(double aMass);
		~DTwoGammaFit_factory(){};
		const string toString(void);
	
	private:
                double fMass;
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

                //DPi0* makePi0(const DLorentzVector& gamma1, const DLorentzVector& gamma2, unsigned int t1, unsigned int t2, oid_t id1, oid_t id2); 
                 //DTwoGammaFit* makeTwoGammaKin(const DKinFit* kfit, const DVector3& vertex);

};

#endif // _DTwoGammaFit_factory_

