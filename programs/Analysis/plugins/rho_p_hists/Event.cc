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
	proton_thrown.SetXYZT(0.0, 0.0, 0.0, 0.0);
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	rho_thrown.m = 0.0;
	rho_thrown.pip.SetXYZT(0.0, 0.0, 0.0, 0.0);
	rho_thrown.pim.SetXYZT(0.0, 0.0, 0.0, 0.0);
	rho_thrown.isfiducial = false;
	rho->Clear(); // delete entries in TClonesArray (without freeing memory)
	Nrho = 0;
}

//------------------
// AddRho
//------------------
void Event::AddRho(TLorentzVector &pip, TLorentzVector &pim)
{
	TClonesArray &rhos = *this->rho;
   rho_t *rho = new(rhos[Nrho++]) rho_t();

	rho->m = (pip+pim).M();
	rho->pip = pip;
	rho->pim = pim;
	rho->isfiducial = IsFiducial(pip) && IsFiducial(pim);
}

//------------------
// IsFiducial
//------------------
bool Event::IsFiducial(TLorentzVector &pion)
{
	double theta = pion.Theta()*TMath::RadToDeg();
	if(theta<2.0 || theta>110.0)return false;
	if(pion.P()<0.500)return false;
	
	return true;
}

