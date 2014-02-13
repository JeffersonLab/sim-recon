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
#include <TDirectoryFile.h>
#include <TLorentzVector.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <TRACKING/DMCThrown.h>
#include <PID/DKinematicData.h>
#include <PID/DChargedTrack.h>
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
	// Create INV_MASS directory
	TDirectory *dir = new TDirectoryFile("INV_MASS","INV_MASS");
	dir->cd();

	// Create histograms
	mm_gp_to_pX = new TH1D("mm_gp_to_pX","Missing mass from #gamma p -> p X",4001, 0.0, 4.0);
	mm_gp_to_pX->SetXTitle("Missing mass (GeV)");

	mass_2gamma = new TH1D("mass_2gamma","2 #gamma invariant mass",4001, 0.0, 4.0);
	mass_2gamma->SetXTitle("Invariant Mass (GeV/c^{2})");

	mass_4gamma = new TH1D("mass_4gamma","4 #gamma invariant mass",4001, 0.0, 4.0);
	mass_4gamma->SetXTitle("Invariant Mass (GeV/c^{2})");

	mass_pip_pim = new TH1D("mass_pip_pim","invariant mass of #pi^{+} and #pi^{-}",4001, 0.0, 4.0);
	mass_pip_pim->SetXTitle("Invariant Mass (GeV/c^{2})");

	t_pX = new TH1D("t_pX","-t for #gammap->pX",20, 0.0, 6.0);
	t_pX->SetXTitle("-t (GeV)");

	sqrt_s = new TH1D("sqrt_s","Center of mass energy #sqrt{s}",2001, 0.0, 6.0);
	sqrt_s->SetXTitle("#sqrt{s}  C.M. energy (GeV)");
	
	// Go back up to the (ROOT) parent directory
	dir->cd("../");
	
	pthread_mutex_init(&mutex, NULL);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_invariant_mass_hists::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects
	vector<const DBeamPhoton*> beam_photons;
	vector<const DPhoton*> photons;
	vector<const DChargedTrack*> chargedtracks;

	loop->Get(beam_photons);	// from truth info
	loop->Get(photons);			// all reconstructed photons (BCAL and FCAL)
	loop->Get(chargedtracks);		// all reconstructed charged (CDC and FDC)

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
	for(unsigned int j=0; j<chargedtracks.size(); j++){
		if(chargedtracks[j]->hypotheses.size()<1)continue;
		const DTrackTimeBased *trk = chargedtracks[j]->hypotheses[0];

		// Use particle mass to decide the type
		int type = 0;
		if(fabs(trk->mass()-0.13957)<0.01){
			type = trk->charge()<0.0 ? 9:8; // initialize to pi-(=9) or pi+(=8)
		}else if(fabs(trk->mass()-0.93827)<0.01 && trk->charge()==1.0){
			type=14;
		}

		// Add TLorentzVector to appropriate container based on charged particle type
		switch(type){
			case 8:	rec_piplus.push_back(MakeTLorentz(trk, 0.13957));		break;
			case 9:	rec_piminus.push_back(MakeTLorentz(trk, 0.13957));	break;
			case 14:	rec_protons.push_back(MakeTLorentz(trk, 0.93827));	break;
			default:
				cout<<"Unknown particle mass: "<<trk->mass()<<" for charge "<<trk->charge()<<endl;
				break;
		}
	} // chargedtracks

	// Some generators don't supply information on the beam photon. If there
	// are no DBeamPhoton objects, then create a 9GeV one here
	if(beam_photons.size()==0){
		DBeamPhoton beam;
		DVector3 mom(0.0, 0.0, 9.0);
		DVector3 pos(0.0, 0.0, 65.0); // center of target
		beam.setMomentum(mom);
		beam.setPosition(pos);
		beam.setMass(0.0);
		beam.setCharge(0.0);
		beam_photons.push_back(&beam);
	}
	
	//--------------------------------------------------------------------
	// Fill histograms below here using values in the rec_XXX containers.

	pthread_mutex_lock(&mutex);
	
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
		for(unsigned int k=0; k<j; k++){
			TLorentzVector sum = rec_photons[j] + rec_photons[k];
			mass_2gamma->Fill(sum.M());
		}
	}

	// 4gamma invariant mass. Loop over all possible combinations
	for(unsigned int j=3; j<rec_photons.size(); j++){
		for(unsigned int k=2; k<j; k++){
			for(unsigned int l=1; l<k; l++){
				for(unsigned int m=0; m<l; m++){
					TLorentzVector sum = rec_photons[j] + rec_photons[k] + rec_photons[l] + rec_photons[m];
					mass_4gamma->Fill(sum.M());
				}
			}
		}
	}

	// pi+, pi- invariant mass. Loop over all possible combinations
	for(unsigned int j=0; j<rec_piplus.size(); j++){
		for(unsigned int k=0; k<rec_piminus.size(); k++){
			mass_pip_pim->Fill( (rec_piplus[j] + rec_piminus[k]).M() );
		}
	}
	
	// Calculate Mandelstam t=(p1-p3)^2
	if(rec_protons.size()==1){
		TLorentzVector beam_photon = beam_photons.size()>0 ? MakeTLorentz(beam_photons[0], 0.0):TLorentzVector(0.0, 0.0, 0.0, 9.0);
		TLorentzVector &proton = rec_protons[0];
		TLorentzVector p3 = beam_photon + target - proton;
		double t = (beam_photon - p3).Mag2();
		t_pX->Fill(-t);
	}

	pthread_mutex_unlock(&mutex);

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
