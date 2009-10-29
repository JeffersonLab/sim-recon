// $Id$
//
//    File: DEventProcessor_phys_tree.cc
// Created: Wed Sep  2 20:25:05 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include "DEventProcessor_phys_tree.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
using namespace std;

#include <TThread.h>
#include <TDirectoryFile.h>
#include <TLorentzVector.h>

#include <TROOT.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include <PID/DKinematicData.h>
#include <PID/DParticle.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>
#include <TRACKING/DMCThrown.h>

#include "Particle.h"

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_phys_tree());
}
} // "C"


//------------------
// DEventProcessor_track_hists
//------------------
DEventProcessor_phys_tree::DEventProcessor_phys_tree()
{
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_track_hists
//------------------
DEventProcessor_phys_tree::~DEventProcessor_phys_tree()
{

}

//------------------
// init
//------------------
jerror_t DEventProcessor_phys_tree::init(void)
{
	// Create PHYSICS directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("PHYSICS");
	if(!dir)dir = new TDirectoryFile("PHYSICS","PHYSICS");
	dir->cd();

	// Here we define a tree with two identical branches based on the Event class.
	// One to hold the thrown values and the other to hold the recon(structed) ones.

	// Create "tree
	tree = new TTree("e","Thrown and Reconstructed Event parameters");

	// Create branches for thrown and reconstructed values
	evt_thrown = new Event();
	evt_recon = new Event();
	tree->Branch("T",&evt_thrown);
	tree->Branch("R",&evt_recon);

	dir->cd("../");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_phys_tree::brun(JEventLoop *loop, int runnumber)
{

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_phys_tree::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects and make TLorentz vectors out of each of them
	vector<const DBeamPhoton*> beam_photons;
	vector<const DMCThrown*> mcthrowns;
	vector<const DPhoton*> photons;
	vector<const DParticle*> particles;

	loop->Get(beam_photons);
	loop->Get(mcthrowns);
	loop->Get(photons);
	loop->Get(particles);

	// Make TLorentzVector for beam photon
	TLorentzVector beam_photon = TLorentzVector(0.0, 0.0, 9.0, 9.0);
	if(beam_photons.size()>0)beam_photon = MakeTLorentz(beam_photons[0], 0.0);
		
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

		// Initially assume it's a pion.
		int type = part->charge()<0.0 ? 9:8; // initialize to pi-(=9) or pi+(=8)

		// If this is a positively charged particle, loop over the thrown protons
		// and see if this matches the momentum to within 10%. If so,
		// set the type of this particle to be a proton
		if(part->charge()>0){
			const TVector3 recon = part->momentum();
			for(unsigned int k=0; k<mcthrowns.size(); k++){
				if(mcthrowns[k]->type!=14)continue; // only interested in thrown protons
				const TVector3 thr = mcthrowns[k]->momentum();
				double dp_over_p = (recon - thr).Mag()/thr.Mag();
				if(fabs(dp_over_p)<0.1){
					type = 14;
					break;
				}
			}
		} 

		// Add TLorentzVector to appropriate container based on charged particle type
		switch(type){
			case 8:	rec_piplus.push_back(MakeTLorentz(part, 0.13957));		break;
			case 9:	rec_piminus.push_back(MakeTLorentz(part, 0.13957));	break;
			case 14:	rec_protons.push_back(MakeTLorentz(part, 0.93827));	break;
		}
	} // particles

	// Create TLorentzVectors for thrown particles
	vector<TLorentzVector> thr_photons;
	vector<TLorentzVector> thr_piplus;
	vector<TLorentzVector> thr_piminus;
	vector<TLorentzVector> thr_protons;
	for(unsigned int k=0; k<mcthrowns.size(); k++){
		switch(mcthrowns[k]->type){
			case  1: thr_photons.push_back(MakeTLorentz(mcthrowns[k], 0.0));		break;
			case  8: thr_piplus.push_back(MakeTLorentz(mcthrowns[k], 0.13957));	break;
			case  9: thr_piminus.push_back(MakeTLorentz(mcthrowns[k], 0.13957));	break;
			case 14: thr_protons.push_back(MakeTLorentz(mcthrowns[k], 0.93827));	break;
		}
	}
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Fill in Event objects for both thrown and reconstructed
	evt_recon->Clear();
	evt_thrown->Clear();
	FillEvent(evt_recon, rec_photons, rec_piplus, rec_piminus, rec_protons);
	FillEvent(evt_thrown, thr_photons, thr_piplus, thr_piminus, thr_protons);
	
	// Copy event number to both trees and add this event to them
	evt_recon->event = eventnumber;
	evt_thrown->event = eventnumber;
	tree->Fill();

	// Unlock mutex
	pthread_mutex_unlock(&mutex);

	return NOERROR;
}


//------------------
// MakeTLorentz
//------------------
TLorentzVector DEventProcessor_phys_tree::MakeTLorentz(const DKinematicData *kd, double mass)
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
// FillEvent
//------------------
void DEventProcessor_phys_tree::FillEvent(Event *evt, vector<TLorentzVector> &photon, vector<TLorentzVector> &pip, vector<TLorentzVector> &pim, vector<TLorentzVector> &proton)
{
	// Add photons
	for(unsigned int i=0; i<photon.size(); i++){
		TClonesArray &prts = *(evt->photon);
		Particle *prt = new(prts[evt->Nphoton++]) Particle();
		prt->p = photon[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add piplus
	for(unsigned int i=0; i<pip.size(); i++){
		TClonesArray &prts = *(evt->pip);
		Particle *prt = new(prts[evt->Npip++]) Particle();
		prt->p = pip[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add piminus
	for(unsigned int i=0; i<pim.size(); i++){
		TClonesArray &prts = *(evt->pim);
		Particle *prt = new(prts[evt->Npim++]) Particle();
		prt->p = pim[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add proton
	for(unsigned int i=0; i<proton.size(); i++){
		TClonesArray &prts = *(evt->proton);
		Particle *prt = new(prts[evt->Nproton++]) Particle();
		prt->p = proton[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}
	
	// Sort particle arrays by energy
	evt->photon->Sort();
	evt->pip->Sort();
	evt->pim->Sort();
	evt->proton->Sort();

	// Calculate W of reconstructed particles
	for(unsigned int i=0; i<photon.size(); i++)evt->W += photon[i];
	for(unsigned int i=0; i<pip.size(); i++)evt->W += pip[i];
	for(unsigned int i=0; i<pim.size(); i++)evt->W += pim[i];
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_phys_tree::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_phys_tree::fini(void)
{

	return NOERROR;
}

