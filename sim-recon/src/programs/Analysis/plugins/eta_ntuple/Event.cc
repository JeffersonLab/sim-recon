// $Id$
//
//    File: Event.cc
// Created: Tue Oct 13 09:55:12 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#include <iostream>
using namespace std;

#include <TVector2.h>

#include "Event.h"

#ifndef _DBG_
#define _DBG_ cerr<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ _DBG_<<endl
#endif

//------------------
// Event (Constructor)
//------------------
Event::Event(){
	fcal = new TClonesArray("fcal_t",100);
	bcal = new TClonesArray("bcal_t",100);
	sc = new TClonesArray("sc_t",100);
}

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
	t = -1000000.0;
	Nstart = 0;
	sc->Clear();	// delete entries in TClonesArray (without freeing memory)
	E_bcal_tot = 0.0;
	Nbcal = 0;
	bcal->Clear();	// delete entries in TClonesArray (without freeing memory)
}

//------------------
// AddFCAL
//------------------
void Event::AddFCAL(TLorentzVector &p, TVector3 &x)
{
	TClonesArray &fcals = *this->fcal;
   fcal_t *fcal = new(fcals[Nfcal++]) fcal_t();

	fcal->p = p;
	fcal->x = x;
}

//------------------
// AddBCAL
//------------------
void Event::AddBCAL(TLorentzVector &p)
{
	TClonesArray &bcals = *this->bcal;
   bcal_t *bcal = new(bcals[Nbcal++]) bcal_t();

	bcal->p = p;
}

//------------------
// AddSC
//------------------
void Event::AddSC(int sector)
{
	TClonesArray &scs = *this->sc;
   sc_t *sc = new(scs[Nstart++]) sc_t();

	sc->phi_center = (float)(sector-1)*2.0*M_PI/40.0;
	
	// Phi angle difference
	float phi_eta = eta_best.Phi();
	TVector2 v1(cos(sc->phi_center), sin(sc->phi_center));
	TVector2 v2(cos(phi_eta), sin(phi_eta));
	sc->phi_diff = v1.DeltaPhi(v2);
	if(sc->phi_diff<0.0)sc->phi_diff+=2.0*M_PI;
}

