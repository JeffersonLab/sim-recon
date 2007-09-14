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
#include "DTwoGammaFit_factory_PI0.h"
#include "JANA/JEvent.h"


//----------------
// Constructor
//----------------
DTwoGammaFit_factory_PI0::DTwoGammaFit_factory_PI0():DTwoGammaFit_factory(0.13498)
{
	// Set defaults

}


//------------------
// evnt
// Pi0 factory: loop over all pair combinations and make pi0
// 		 regardless of the pair-mass at this point (MK)
//------------------
/*jerror_t DTwoGammaFit_factory::evnt(JEventLoop *eventLoop, int eventnumber)
{
	vector<const DPhoton*> photons;
	eventLoop->Get(photons);
	

	// Loop over all photons and fit candidates 
        for (unsigned int i = 0; i < photons.size() ; i++) {
//           if (photons[i]->getTag() == 2) continue; // apply charged clusters cut here
          for (unsigned int j = i+1; j < photons.size() ; j++) {
//              if (photons[j]->getTag() == 2) continue; // apply charged clusters cut here
                
                DKinFit *kfit = new DKinFit(); 
                vector<const DKinematicData*> kindata;
                kfit->SetVerbose(0);

                kindata.clear();
                kindata.push_back((DKinematicData*)photons[i]);
                kindata.push_back((DKinematicData*)photons[j]);

                kfit->SetFinal(kindata);
//                kfit->FitTwoGammas(0.13498); 
                cout << " Fitting two photons to:  "  << fMass << endl;
                kfit->FitTwoGammas(fMass); 
 
                vector<DKinematicData*> kinout = kfit->GetFinal_out();

                DTwoGammaFit* fit2g = new DTwoGammaFit();

 get final pair kinematics 
                const DLorentzVector& P41 = kfit->FitP4(0);
                const DLorentzVector& P42 = kfit->FitP4(1);
                DLorentzVector p4 = P41 + P42; 
               
                if (kinout.size() != 2) continue;
                TLorentzVector P4 = kinout[0]->lorentzMomentum() + kinout[1]->lorentzMomentum(); 
		DVector3 vertex = 0.5*(photons[i]->position() + photons[j]->position());

// set two gamma kinematics
                fit2g->setMass( P4.M() );
                fit2g->setMomentum( P4.Vect() );
                fit2g->setPosition( vertex );
                fit2g->setProb( kfit->Prob());
                fit2g->setChi2( kfit->Chi2());
                fit2g->setPulls( kfit->GetPull(1),  1);
                //fit2g->setChildFits(??,0);

		_data.push_back( fit2g );

           } 
        } 


	return NOERROR;
}
*/
//------------------
// toString
//------------------
const string DTwoGammaFit_factory_PI0::toString(void)
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
//		printcol("%5.2f", pi0s->getM());
//		printcol("%4.0f", fcalhit->t);
		printrow();
	}

	return _table;
}

