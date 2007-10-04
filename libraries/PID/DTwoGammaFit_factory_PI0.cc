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

