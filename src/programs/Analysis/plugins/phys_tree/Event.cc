// $Id$
//
//    File: Event.cc
// Created: Tue Oct 13 09:55:12 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//



#include "Event.h"
#include "Particle.h"

//------------------
// Event (ctor)
//------------------
Event::Event()
{
	pip = new TClonesArray("Particle", MAX_PART);
	pim = new TClonesArray("Particle", MAX_PART);
	proton = new TClonesArray("Particle", MAX_PART);
	photon = new TClonesArray("Particle", MAX_PART);
}

//------------------
// ~Event (dtor)
//------------------
Event::~Event()
{
	delete pip;
	delete pim;
	delete proton;
	delete photon;
}

//------------------
// Clear
//------------------
void Event::Clear(void)
{
	Npip = 0;
	Npim = 0;
	Nproton = 0;
	Nphoton = 0;
	pip->Clear();	// delete entries in TClonesArray (without freeing memory)
	pim->Clear();	// delete entries in TClonesArray (without freeing memory)
	proton->Clear();	// delete entries in TClonesArray (without freeing memory)
	photon->Clear();	// delete entries in TClonesArray (without freeing memory)
	
	target.SetXYZT(0.0, 0.0, 0.0, 0.0);
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	W.SetXYZT(0.0, 0.0, 0.0, 0.0);

	all_mesons_fiducial = false;
	all_photons_fiducial = false;
	all_protons_fiducial = false;
}

#if 0
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

#endif
