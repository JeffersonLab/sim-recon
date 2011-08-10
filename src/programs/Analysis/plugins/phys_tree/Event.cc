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
	Kp = new TClonesArray("Particle", MAX_PART);
	Km = new TClonesArray("Particle", MAX_PART);
	proton = new TClonesArray("Particle", MAX_PART);
	photon = new TClonesArray("Particle", MAX_PART);
	neutron = new TClonesArray("Particle", MAX_PART);

	pip_match = new TClonesArray("Particle", MAX_PART);
	pim_match = new TClonesArray("Particle", MAX_PART);
	Kp_match = new TClonesArray("Particle", MAX_PART);
	Km_match = new TClonesArray("Particle", MAX_PART);
	proton_match = new TClonesArray("Particle", MAX_PART);
	photon_match = new TClonesArray("Particle", MAX_PART);
	neutron_match = new TClonesArray("Particle", MAX_PART);
}

//------------------
// ~Event (dtor)
//------------------
Event::~Event()
{
	delete pip;
	delete pim;
	delete Kp;
	delete Km;
	delete proton;
	delete photon;
	delete neutron;

	delete pip_match;
	delete pim_match;
	delete Kp_match;
	delete Km_match;
	delete proton_match;
	delete photon_match;
	delete neutron_match;
}

//------------------
// Clear
//------------------
void Event::Clear(void)
{
	Npip = 0;
	Npim = 0;
	NKp = 0;
	NKm = 0;
	Nproton = 0;
	Nphoton = 0;
	Nneutron = 0;
	pip->Clear();		// delete entries in TClonesArray (without freeing memory)
	pim->Clear();		// delete entries in TClonesArray (without freeing memory)
	Kp->Clear();		// delete entries in TClonesArray (without freeing memory)
	Km->Clear();		// delete entries in TClonesArray (without freeing memory)
	proton->Clear();	// delete entries in TClonesArray (without freeing memory)
	photon->Clear();	// delete entries in TClonesArray (without freeing memory)
	neutron->Clear();	// delete entries in TClonesArray (without freeing memory)

	pip_match->Clear();
	pim_match->Clear();
	Kp_match->Clear();
	Km_match->Clear();
	proton_match->Clear();
	photon_match->Clear();
	neutron_match->Clear();
	
	target.SetXYZT(0.0, 0.0, 0.0, 0.0);
	beam.SetXYZT(0.0, 0.0, 0.0, 0.0);
	vertex.SetXYZ(0.0, 0.0, 65.0);
	W.SetXYZT(0.0, 0.0, 0.0, 0.0);

	all_fiducial = false;
	all_mesons_fiducial = false;
	all_photons_fiducial = false;
	all_neutrons_fiducial = false;
	all_protons_fiducial = false;
}


