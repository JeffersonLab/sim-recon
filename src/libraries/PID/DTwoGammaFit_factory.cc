// $Id: DTwoGammaFit_factory.cc 2496 2007-03-12 00:35:46Z kornicer $
//
//    File: DTwoGammaFit_factory.cc
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#include <math.h>

#include <JANA/JEvent.h>
using namespace jana;

#include "DTwoGammaFit.h"
#include "DKinFit.h"
#include "DPhoton.h"
#include "DTwoGammaFit_factory.h"


bool SortByProb(const DTwoGammaFit *a,const DTwoGammaFit *b){
  return (a->getProb() > b->getProb());
}


//----------------
// Constructor
//----------------
DTwoGammaFit_factory::DTwoGammaFit_factory(double aMass)
{
	// Set defaults
        if (aMass > 0) { 
           fMass = aMass; 
         }
         else {
           cout << "Mass undefined, set to pi0" << endl; 
           fMass = 0.135; 
         }

}


//------------------
// evnt
// TwoGammaGit factory: loop over all pair combinations and fit them to the supplied mass using
//			kinemati fitter. Tags, like PI0, ETA (and ETAP eventually) can be used in derived 
//			factories.
//------------------
jerror_t DTwoGammaFit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DPhoton*> photons;
	eventLoop->Get(photons);
   
        JObject::oid_t nPairs=0;
	// Loop over all fit candidates.
        for (unsigned int i = 0; i < photons.size() ; i++) {
            if ( photons[i]->getTag() == DPhoton::kCharge  ) continue;
            for (unsigned int j = i+1; j < photons.size() ; j++) {
                if ( photons[j]->getTag() == DPhoton::kCharge  ) continue;

                nPairs++;
                DKinFit *kfit = new DKinFit(); 
                kfit->SetVerbose(0);

                vector<DKinematicData> kindata;
                kindata.clear();
                kindata.push_back( *photons[i] );
                kindata.push_back( *photons[j] );

                kfit->SetFinal(kindata);
                kfit->FitTwoGammas(fMass,1.); // second parameter is scale for photon error matrix
 
                vector<DKinematicData> kinout = kfit->GetFinal_out();
                if (kinout.size() != 2) continue;

                DTwoGammaFit* fit2g = new DTwoGammaFit( nPairs );

// set two gamma kinematics
                DLorentzVector P4U = photons[i]->lorentzMomentum() + photons[j]->lorentzMomentum(); 
                DLorentzVector P4 = kinout[0].lorentzMomentum() + kinout[1].lorentzMomentum(); 
		DVector3 vertex = 0.5*(photons[i]->position() + photons[j]->position());

                fit2g->setUMass( P4U.M() );
                fit2g->setMass( P4.M() );
                fit2g->setMomentum( P4.Vect() );
                fit2g->setPosition( vertex );
                fit2g->setProb( kfit->Prob());
                fit2g->setChi2( kfit->Chi2());
					 fit2g->setNdf( kfit->Ndf());
                for (int k=0; k < 6; k++) {
                   fit2g->setPulls( kfit->GetPull(k),  k);
                }
                fit2g->setChildTag( photons[i]->getTag(), 0);
                fit2g->setChildTag( photons[j]->getTag(), 1);
                fit2g->setChildID( photons[i]->id, 0);
                fit2g->setChildID( photons[j]->id, 1);
                fit2g->setChildMom(  photons[i]->lorentzMomentum() , 0);
                fit2g->setChildMom(  photons[j]->lorentzMomentum() , 1);

                fit2g->setChildFit(  kinout[0] , 0);
                fit2g->setChildFit(  kinout[1] , 1);
 

		_data.push_back( fit2g );
 
		delete kfit;

           } 
        } 

	// sort by probability so most probable fits are listed first
	sort(_data.begin(), _data.end(), SortByProb);

	return NOERROR;
}


