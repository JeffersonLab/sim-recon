// $Id: DPi0_factory.cc 2496 2007-03-12 00:35:46Z kornicer $
//
//    File: DPi0_factory.cc
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#include <math.h>

#include "DPi0.h"
#include "DPhoton.h"
#include "DPi0_factory.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DPi0_factory::DPi0_factory()
{
	// Set defaults

}


//------------------
// evnt
// Pi0 factory: loop over all pair combinations and make pi0
// 		 regardless of the pair-mass at this point (MK)
//------------------
jerror_t DPi0_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DPhoton*> photons;
	eventLoop->Get(photons);
	
	// Loop over all photons and make pi0 candidates 
        for (unsigned int i = 0; i < photons.size() ; i++) {

            if (photons[i]->getdThetaCharge() < 0.05 ) continue;
            for (unsigned int j = i+1; j < photons.size() ; j++) {

                if (photons[j]->getdThetaCharge() < 0.05 ) continue;
                
		DPi0 *pi0 =  makePi0( photons[i], photons[j] );

		_data.push_back(pi0);

           } 
        } 

	return NOERROR;
}


// create pi0 candidate from two photons 
DPi0* DPi0_factory::makePi0(const DPhoton* gamma1, const DPhoton* gamma2) 
{

        DPi0* pi0 = new DPi0;
        DLorentzVector mom4 = gamma1->lorentzMomentum() + gamma2->lorentzMomentum();
        DVector3 vertex = 0.5*(gamma1->position() + gamma2->position());

        pi0->setMass( mom4.M() );
        pi0->setMomentum( mom4.Vect() );
        pi0->setPosition( vertex );
        pi0->setChildrenTag(gamma1->getTag(), gamma1->getTag());
        pi0->setChildrenID(gamma1->id, gamma2->id);

        return pi0;
}


