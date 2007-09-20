// $Id: DTwoGammaFit_factory.cc 2496 2007-03-12 00:35:46Z kornicer $
//
//    File: DTwoGammaFit_factory.cc
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#include <math.h>

#include "DTwoGammaFit.h"
#include "DKinFit.h"
#include "DPhoton.h"
#include "DTwoGammaFit_factory.h"
#include "JANA/JEvent.h"


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
   
	// Loop over all fit candidates.
        for (unsigned int i = 0; i < photons.size() ; i++) {
            if ( photons[i]->getdThetaCharge() < 0.05  ) continue;
            for (unsigned int j = i+1; j < photons.size() ; j++) {
                if ( photons[j]->getdThetaCharge() < 0.05  ) continue;

                
                DKinFit *kfit = new DKinFit(); 
                vector<const DKinematicData*> kindata;
                kfit->SetVerbose(0);

                kindata.clear();
                kindata.push_back((DKinematicData*)photons[i]);
                kindata.push_back((DKinematicData*)photons[j]);

                kfit->SetFinal(kindata);
                kfit->FitTwoGammas(fMass,1.); // second parameter is scale for photon error matrix
 
                vector<DKinematicData*> kinout = kfit->GetFinal_out();
                if (kinout.size() != 2) continue;

                DTwoGammaFit* fit2g = new DTwoGammaFit();

                //TLorentzVector P4 = kinout[0]->lorentzMomentum() + kinout[1]->lorentzMomentum(); 
                TLorentzVector P4 = photons[i]->lorentzMomentum() + photons[j]->lorentzMomentum(); 
		DVector3 vertex = 0.5*(photons[i]->position() + photons[j]->position());

// set two gamma kinematics
                fit2g->setMass( P4.M() );
                fit2g->setMomentum( P4.Vect() );
                fit2g->setPosition( vertex );
                fit2g->setProb( kfit->Prob());
                fit2g->setChi2( kfit->Chi2());
                for (int k=0; k < 6; k++) {
                   fit2g->setPulls( kfit->GetPull(k),  k);
                }
                fit2g->setChildTag( (int) photons[i]->getTag(), 0);
                fit2g->setChildTag( (int) photons[j]->getTag(), 1);
                fit2g->setChildID( photons[i]->id, 0);
                fit2g->setChildID( photons[j]->id, 1);
                fit2g->setChildFit(  *kinout[0] , 0);
                fit2g->setChildFit(  *kinout[1] , 1);
 
		_data.push_back( fit2g );

           } 
        } 


	return NOERROR;
}

//------------------
// toString
//------------------
const string DTwoGammaFit_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	Get();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row:   E(GeV):   Px(cm):   Py(cm):   Pz(cm):   M(GeV):   Chi2:");
	
	for(unsigned int i=0; i<_data.size(); i++){
		DTwoGammaFit *pions = _data[i];

		printnewrow();
		printcol("%d",	i);
		printcol("%6.2f", pions->energy());
		printcol("%6.2f", pions->lorentzMomentum().X());
		printcol("%6.2f", pions->lorentzMomentum().Y());
		printcol("%6.2f", pions->lorentzMomentum().Z());
		printcol("%6.2f", pions->mass());
		printcol("%6.2f", pions->getChi2());
		printrow();
	}

	return _table;
}

