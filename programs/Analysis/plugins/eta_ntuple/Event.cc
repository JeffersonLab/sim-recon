// $Id$
//
//    File: Event.cc
// Created: Tue Oct 13 09:55:12 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//



#include "Event.h"


//------------------
// Clear
//------------------
void Event::Clear(void)
{
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	proton_thrown.SetXYZT(0.0, 0.0, 0.0, 0.0);
	eta_thrown.SetXYZT(0.0, 0.0, 0.0, 0.0);
	eta_best.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	prod_mech = 0;
	decay_mode = 0;
	fcal->Clear(); // delete entries in TClonesArray (without freeing memory)
	Nfcal = 0;
}

//------------------
// AddRho
//------------------
void Event::AddFCAL(TLorentzVector &p, TVector3 &x)
{
	TClonesArray &fcals = *this->fcal;
   fcal_t *fcal = new(fcals[Nfcal++]) fcal_t();

	fcal->p = p;
	fcal->x = x;
}

