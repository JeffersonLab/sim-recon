// $Id: DPIDPi0_factory.cc 2496 2007-03-12 00:35:46Z kornicer $
//
//    File: DPIDPi0_factory.cc
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#include <math.h>

#include "DPIDPi0.h"
#include "DPIDPhoton.h"
#include "DPIDPi0_factory.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DPIDPi0_factory::DPIDPi0_factory()
{
	// Set defaults

}


//------------------
// evnt
// Pi0 facttory: loop over all pair combinations and make pi0
// 		 regatdless of the pair-mass at this point (MK)
//------------------
jerror_t DPIDPi0_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DPIDPhoton*> photons;
	eventLoop->Get(photons);
	
	// Loop over all photons and make pi0 candidates 
        for (unsigned int i = 0; i < photons.size() ; i++) {
          for (unsigned int j = i+1; j < photons.size() ; j++) {

                TLorentzVector P1 = photons[i]->getMom4();
                TLorentzVector P2 = photons[j]->getMom4();
                const unsigned int tag1 = photons[i]->getTag();
                const unsigned int tag2 = photons[j]->getTag();
		DPIDPi0 *pi0 =  makePi0(P1, P2, tag1, tag2);

		_data.push_back(pi0);

           } 
        } 

	return NOERROR;
}

//------------------
// toString
//------------------
const string DPIDPi0_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   E(GeV):   x(cm):   y(cm):   z(cm):	M(GeV):   t(ns):");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DPIDPi0 *pions = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%5.2f", pions->getMom4().T());
		printcol("%5.2f", pions->getMom4().X());
		printcol("%5.2f", pions->getMom4().Y());
		printcol("%5.2f", pions->getMom4().Z());
//		printcol("%5.2f", pi0s->getM());
//		printcol("%4.0f", fcalhit->t);
		printrow();
	}

	return _table;
}

// create pi0 candidate from two photons 
DPIDPi0* DPIDPi0_factory::makePi0(const TLorentzVector gamma1, const TLorentzVector gamma2, const unsigned int tag1, const unsigned int tag2) 
{

        DPIDPi0* pi0 = new DPIDPi0;
        
        pi0->setMom4(gamma1, gamma2);
        pi0->setOrig(tag1, tag2);

        return pi0;
}

