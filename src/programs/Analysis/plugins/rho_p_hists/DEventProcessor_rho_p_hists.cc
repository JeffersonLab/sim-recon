// $Id: DEventProcessor_rho_p_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_rho_p_hists.cc
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
#include <PID/DParticle.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>

#include "DEventProcessor_rho_p_hists.h"

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
	app->AddProcessor(new DEventProcessor_rho_p_hists());
}
}


//------------------
// init
//------------------
jerror_t DEventProcessor_rho_p_hists::init(void)
{
	// Create histograms
	mm_gp_to_pX = new TH1D("mm_gp_to_pX","Missing mass from #gamma p -> p X",4001, 0.0, 4.0);
	mm_gp_to_pX->SetXTitle("Missing mass (GeV)");
	mm_gp_to_pX_thrown = (TH1D*)mm_gp_to_pX->Clone("mm_gp_to_pX_thrown");

	t_pX = new TH1D("t_pX","-t for #gammap->pX",20, 0.0, 6.0);
	t_pX->SetXTitle("-t (GeV)");

	sqrt_s = new TH1D("sqrt_s","Center of mass energy #sqrt{s}",2001, 0.0, 6.0);
	sqrt_s->SetXTitle("#sqrt{s}  C.M. energy (GeV)");
	
	// Create tree
	evt = new Event();
	tree = new TTree("t", "#rho p events");
	tree->Branch("E",&evt);
	
	pthread_mutex_init(&mutex, NULL);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_rho_p_hists::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects
	vector<const DBeamPhoton*> beam_photons;
	vector<const DParticle*> particles;
	vector<const DParticle*> particles_thrown;

	loop->Get(beam_photons);	// from truth info
	loop->Get(particles);		// all reconstructed charged (CDC and FDC)
	loop->Get(particles_thrown, "THROWN");		// all thrown charged (CDC and FDC)

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Create TLorentzVectors for reconstructed charged particles
	vector<TLorentzVector> rec_piplus;
	vector<TLorentzVector> rec_piminus;
	vector<TLorentzVector> rec_protons;
	SortChargedParticles(particles, rec_piplus, rec_piminus, rec_protons);

	// Create TLorentzVectors for thrown charged particles
	vector<TLorentzVector> thrown_piplus;
	vector<TLorentzVector> thrown_piminus;
	vector<TLorentzVector> thrown_protons;
	SortChargedParticles(particles_thrown, thrown_piplus, thrown_piminus, thrown_protons);

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
	
	evt->Clear();
	evt->event = eventnumber;
	
	// Loop over beam photons and fill histos for each "tagged" photon for this event
	for(unsigned int i=0; i<beam_photons.size(); i++){

		// Make TLorentzVector for beam photon
		TLorentzVector beam_photon = MakeTLorentz(beam_photons[i], 0.0);
		evt->beam = beam_photon;

		// Center of mass energy
		sqrt_s->Fill((beam_photon+target).M()); // Fill Center of Mass energy histo
				
		// Missing mass from gamma + p -> p + X
		if(rec_protons.size()==1){
			TLorentzVector missing = beam_photon + target - rec_protons[0];
			mm_gp_to_pX->Fill(missing.M());
		}

		// Missing mass from gamma + p -> p + X
		if(thrown_protons.size()==1){
			TLorentzVector missing = beam_photon + target - thrown_protons[0];
			mm_gp_to_pX_thrown->Fill(missing.M());
		}

	} // beam_photons


	//-------------------------------------------------------------------
	// Below here are the histos that do not depend on the tagged photon
	// energy and so should be filled outside of the above loop.

	// We are only interested in single rho events so return now if there
	// are not exactly one thrown pi+ and one thrown pi-
	if(thrown_piplus.size()!=1)return NOERROR;
	if(thrown_piminus.size()!=1)return NOERROR;

	// Determine whether both thrown pions are fiducial
	evt->rho_thrown.isfiducial = IsFiducial(thrown_piplus[0]) && IsFiducial(thrown_piminus[0]);
	
	// Thrown pion parameters
	evt->rho_thrown.pip = thrown_piplus[0];
	evt->rho_thrown.pim = thrown_piminus[0];
	evt->rho_thrown.m = (evt->rho_thrown.pip+evt->rho_thrown.pim).M();

	// pi+, pi- invariant mass. Loop over all possible combinations
	for(unsigned int j=0; j<rec_piplus.size(); j++){
		for(unsigned int k=0; k<rec_piminus.size(); k++){
			evt->AddRho(rec_piplus[j], rec_piminus[k]);
		}
	}
	
	// Thrown proton
	if(thrown_protons.size()==1)evt->proton_thrown = thrown_protons[0];

	// Calculate Mandelstam t=(p1-p3)^2
	if(rec_protons.size()==1){
		TLorentzVector beam_photon = beam_photons.size()>0 ? MakeTLorentz(beam_photons[0], 0.0):TLorentzVector(0.0, 0.0, 0.0, 9.0);
		TLorentzVector &proton = rec_protons[0];
		TLorentzVector p3 = beam_photon + target - proton;
		double t = (beam_photon - p3).Mag2();
		t_pX->Fill(-t);
	}
	
	// Fill tree
	evt->rho->Sort(); // sort by closeness of invariant mass to 770MeV/c^2 (uses rho_t::Compare );
	tree->Fill();

	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// SortChargedParticles
//------------------
void DEventProcessor_rho_p_hists::SortChargedParticles(vector<const DParticle*> &particles, vector<TLorentzVector> &rec_piplus, vector<TLorentzVector> &rec_piminus, vector<TLorentzVector> &rec_protons)
{

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
}

//------------------
// MakeTLorentz
//------------------
TLorentzVector DEventProcessor_rho_p_hists::MakeTLorentz(const DKinematicData *kd, double mass)
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
// IsFiducial
//------------------
bool DEventProcessor_rho_p_hists::IsFiducial(TLorentzVector &pion)
{
	double theta = pion.Theta()*TMath::RadToDeg();
	if(theta<2.0 || theta>110.0)return false;
	if(pion.P()<0.500)return false;
	
	return true;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_rho_p_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_rho_p_hists::fini(void)
{
	// Histograms are automatically written to file by hd_root program (or equivalent)
	// so we don't need to do anything here.

	return NOERROR;
}
