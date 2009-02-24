// $Id: DEventProcessor_invariant_mass_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_invariant_mass_hists.cc
// Created: Thur Jan 11 15:42:21 EDT 2005
// Creator: davidl 
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>
#include <TLorentzVector.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <PID/DKinematicData.h>
#include <PID/DParticle.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>

#include "DEventProcessor_invariant_mass_hists.h"

//==========================================================================
// This file contains code for a plugin that will create a few histograms
// based on the reconstructed data. It can mainly serve as an example
// of how one could access the reconstructed data in their own analysis
// program.
//
// Histograms are created in the init() method and they are filled in
// the evnt() method.
//==========================================================================

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_invariant_mass_hists());
}
}


//------------------
// init
//------------------
jerror_t DEventProcessor_invariant_mass_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();

	// Create THROWN directory
	TDirectory *dir = new TDirectoryFile("INV_MASS","INV_MASS");
	dir->cd();

	// Create histograms
	mm_gp_to_pX = new TH1D("mm_gp_to_pX","Missing mass from #gamma p -> p X",4001, 0.0, 4.0);
	mm_gp_to_pX->SetXTitle("Missing mass (GeV)");

	mass_2gamma = new TH1D("mass_2gamma","2 #gamma invariant mass",4001, 0.0, 4.0);
	mass_2gamma->SetXTitle("Invariant Mass (GeV/c^{2})");

	mass_pip_pim = new TH1D("mass_pip_pim","invariant mass of #pi^{+} and #pi^{-}",4001, 0.0, 4.0);
	mass_pip_pim->SetXTitle("Invariant Mass (GeV/c^{2})");

	sqrt_s = new TH1D("sqrt_s","Center of mass energy #sqrt{s}",2001, 0.0, 6.0);
	sqrt_s->SetXTitle("#sqrt{s}  C.M. energy (GeV)");
	
	// Go back up to the (ROOT) parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_invariant_mass_hists::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects and make TLorentz vectors out of each of them
	vector<const DBeamPhoton*> beam_photons;
	vector<const DPhoton*> photons;
	vector<const DParticle*> particles;

	loop->Get(beam_photons);
	loop->Get(photons);
	loop->Get(particles);

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Create TLorentzVectors for reconstructed photons
	vector<TLorentzVector> rec_photons;
	for(unsigned int j=0; j<photons.size(); j++)rec_photons.push_back(MakeTLorentz(photons[j], 0.0));	

	// Create TLorentzVectors for reconstructed charged particles
	vector<TLorentzVector> rec_piplus;
	vector<TLorentzVector> rec_piminus;
	vector<TLorentzVector> rec_protons;

	// Loop over charged particles turning them into TLorentzVector objects
	// and sorting them into various containers declared just above.
	for(unsigned int j=0; j<particles.size(); j++){
		const DParticle *part = particles[j];

		// If this came from the HDParSim factory it will have a DMCThrown object
		// associated with it that we can use to get the particle type since
		// we currently don't have that info in reconstruction. If it is not there,
		// then assume it's a pion.
		int type = part->charge()<0.0 ? 9:8; // initialize to pi-(=9) or pi+(=8)

		// Here we try and get the "truth" object DMCThrown. The access mechanism
		// forces us to get it as a list, but there should be at most 1.
		vector<const DMCThrown*> throwns;
		part->Get(throwns);
		// if DMCThrown was found, overwrite pion with actual particle type
		if(throwns.size()>0)type = throwns[0]->type;			

		// Add TLorentzVector to appropriate container based on charged particle type
		switch(type){
			case 8:	rec_piplus.push_back(MakeTLorentz(part, 0.13957));		break;
			case 9:	rec_piminus.push_back(MakeTLorentz(part, 0.13957));	break;
			case 14:	rec_protons.push_back(MakeTLorentz(part, 0.93827));	break;
		}
	} // particles

	
	//--------------------------------------------------------------------
	// Fill histograms below here using values in the rec_XXX constainers.


	// Loop over beam photons and fill histos for each "tagged" photon for this event
	for(unsigned int i=0; i<beam_photons.size(); i++){

		// Make TLorentzVector for beam photon
		TLorentzVector beam_photon = MakeTLorentz(beam_photons[i], 0.0);

		// Center of mass energy
		sqrt_s->Fill((beam_photon+target).M()); // Fill Center of Mass energy histo
				
		// Missing mass from gamma + p -> p + X
		if(rec_protons.size()==1){
			TLorentzVector missing = beam_photon + target - rec_protons[0];
			mm_gp_to_pX->Fill(missing.M());
		}

	} // beam_photons


	//-------------------------------------------------------------------
	// Below here are the histos that do not depend on the tagged photon
	// energy and so should be filled outside of the above loop.

	// 2gamma invariant mass. Loop over all possible combinations
	for(unsigned int j=1; j<rec_photons.size(); j++){
		TLorentzVector &ph1 = rec_photons[j];
		for(unsigned int k=0; k<j; k++){
			TLorentzVector &ph2 = rec_photons[k];
			
			TLorentzVector sum = ph1 + ph2;
			mass_2gamma->Fill(sum.M());
		}
	}

	// pi+, pi- invariant mass. Loop over all possible combinations
	for(unsigned int j=0; j<rec_piplus.size(); j++){
		for(unsigned int k=0; k<rec_piminus.size(); k++){
			mass_pip_pim->Fill( (rec_piplus[j] + rec_piminus[k]).M() ); // (a more compact way than the 2 gamma example above)
		}
	}		

	return NOERROR;
}

//------------------
// MakeTLorentz
//------------------
TLorentzVector DEventProcessor_invariant_mass_hists::MakeTLorentz(const DKinematicData *kd, double mass)
{
	// Create a ROOT TLorentzVector object out of a Hall-D DKinematic Data object.
	// Here, we have the mass passed in explicitly rather than use the mass contained in
	// the DKinematicData object itself. This is because right now (Feb. 2009) the
	// PID code is not mature enough to give reasonable guesses. See above code.

	double p = kd->momentum().Mag();
	double theta = kd->momentum().Theta();
	double phi = kd->momentum().Phi();
	double px = p*sin(theta)*cos(phi);
	double py = p*sin(theta)*sin(phi);
	double pz = p*cos(theta);
	double E = sqrt(mass*mass + p*p);
	
	return TLorentzVector(px,py,pz,E);
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_invariant_mass_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_invariant_mass_hists::fini(void)
{
	// Histograms are automatically written to file by hd_root program (or equivalent)
	// so we don't need to do anything here.

	return NOERROR;
}
