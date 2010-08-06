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
#include <PID/DChargedTrack.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>
#include <PID/DPhysicsEvent.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackTimeBased.h>
#include "Particle.h"

bool static CompareLorentzEnergy(const TLorentzVector &a, const TLorentzVector &b){
  return a.E()<b.E();
}

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
	tree_thrwn = new TTree("thrown","Thrown Event parameters");
	tree_recon = new TTree("recon","Reconstructed Event parameters");

	// Create branches for thrown and reconstructed values
	evt_thrown = new Event();
	evt_recon = new Event();
	tree_thrwn->Branch("T",&evt_thrown);
	tree_recon->Branch("R",&evt_recon);
	
	// Empty tree with thrown and recon values as friends
	TTree *evt = new TTree("evt","Thrown and Reconstructed Event parameters");
	evt->AddFriend("tree_thrwn");
	evt->AddFriend("tree_recon");

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
	vector<const DChargedTrack*> chargedtracks;
	vector<const DPhysicsEvent*> physicsevents;

	loop->Get(beam_photons);
	loop->Get(mcthrowns);
	loop->Get(photons);
	loop->Get(chargedtracks);
	loop->Get(physicsevents);

	// Make TLorentzVector for beam photon
	TLorentzVector beam_photon = TLorentzVector(0.0, 0.0, 9.0, 9.0);
	if(beam_photons.size()>0)beam_photon = MakeTLorentz(beam_photons[0], 0.0);
		
	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Find the DPhysicsEvent with the most particles and only use that one.
	// This is not a long term solution, but is motivated by the fact that
	// we have only one set of DMCThrown particles and one DBeamPhoton
	const DPhysicsEvent *physicsevent = NULL;
	int max_parts=0;
	for(unsigned int i=0; i<physicsevents.size(); i++){
		const DPhysicsEvent *pe = physicsevents[i];
		int Nparts = pe->pip.size() + pe->pim.size()
		        + pe->photon.size() + pe->proton.size()
		        + pe->Kp.size()     + pe->Km.size()
		        + pe->otherp.size() + pe->otherm.size();
		if(Nparts>max_parts || physicsevent==NULL){
			physicsevent = pe;
			max_parts = Nparts;
		}
	}
	
	// Created TLorentzVector objects for each of the common particle types
	particle_set rec;
	for(unsigned int j=0; j<physicsevent->photon.size(); j++){
		// photon
		rec.photons.push_back(MakeTLorentz(physicsevent->photon[j], 0.0));
	}
	for(unsigned int j=0; j<physicsevent->pip.size(); j++){
		// pi+
		rec.piplus.push_back(MakeTLorentz(physicsevent->pip[j], 0.13957));
	}
	for(unsigned int j=0; j<physicsevent->pim.size(); j++){
		// pi-
		rec.piminus.push_back(MakeTLorentz(physicsevent->pim[j], 0.13957));
	}
	for(unsigned int j=0; j<physicsevent->proton.size(); j++){
		// proton
		rec.protons.push_back(MakeTLorentz(physicsevent->proton[j], 0.93827));
	}
	for(unsigned int j=0; j<physicsevent->Kp.size(); j++){
		// K+
		rec.Kplus.push_back(MakeTLorentz(physicsevent->Kp[j], 0.493677));
	}
	for(unsigned int j=0; j<physicsevent->Km.size(); j++){
		// Ki+
		rec.Kminus.push_back(MakeTLorentz(physicsevent->Km[j], 0.493677));
	}

	// Create TLorentzVectors for thrown particles
	particle_set thr;
	for(unsigned int k=0; k<mcthrowns.size(); k++){
		switch(mcthrowns[k]->type){
			case  1: thr.photons.push_back(MakeTLorentz(mcthrowns[k], 0.0));		break;
			case  8: thr.piplus.push_back(MakeTLorentz(mcthrowns[k], 0.13957));		break;
			case  9: thr.piminus.push_back(MakeTLorentz(mcthrowns[k], 0.13957));	break;
			case 11: thr.Kplus.push_back(MakeTLorentz(mcthrowns[k], 0.493677));		break;
			case 12: thr.Kminus.push_back(MakeTLorentz(mcthrowns[k], 0.493677));	break;
			case 14: thr.protons.push_back(MakeTLorentz(mcthrowns[k], 0.93827));	break;
		}
	}
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Fill in Event objects for both thrown and reconstructed
	evt_recon->Clear();
	evt_thrown->Clear();
	FillEvent(evt_recon, rec, thr);
	FillEvent(evt_thrown, thr, rec);
	
	// Copy event number to both trees and add this event to them
	evt_recon->event = eventnumber;
	evt_thrown->event = eventnumber;
	tree_thrwn->Fill();
	tree_recon->Fill();

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
void DEventProcessor_phys_tree::FillEvent(Event *evt, particle_set &pset, particle_set &pset_match)
{
	vector<TLorentzVector> &photon = pset.photons;
	vector<TLorentzVector> &pip = pset.piplus;
	vector<TLorentzVector> &pim = pset.piminus;
	vector<TLorentzVector> &proton = pset.protons;

	vector<TLorentzVector> &photon_match = pset_match.photons;
	vector<TLorentzVector> &pip_match = pset_match.piplus;
	vector<TLorentzVector> &pim_match = pset_match.piminus;
	vector<TLorentzVector> &proton_match = pset_match.protons;

	// Sort particle arrays by energy
	sort(photon.begin(), photon.end(), CompareLorentzEnergy);
	sort(pip.begin(), pip.end(), CompareLorentzEnergy);
	sort(pim.begin(), pim.end(), CompareLorentzEnergy);
	sort(proton.begin(), proton.end(), CompareLorentzEnergy);

	// Add photons
	for(unsigned int i=0; i<photon.size(); i++){
		TClonesArray &prts_match = *(evt->photon_match);
		Particle *prt_match = new(prts_match[evt->Nphoton]) Particle();
		prt_match->p = FindBestMatch(photon[i], photon_match);
		prt_match->x.SetXYZ(0,0,65); // FIXME!!!

		TClonesArray &prts = *(evt->photon);
		Particle *prt = new(prts[evt->Nphoton++]) Particle();
		prt->p = photon[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add piplus
	for(unsigned int i=0; i<pip.size(); i++){
		TClonesArray &prts_match = *(evt->pip_match);
		Particle *prt_match = new(prts_match[evt->Npip]) Particle();
		prt_match->p = FindBestMatch(pip[i], pip_match);
		prt_match->x.SetXYZ(0,0,65); // FIXME!!!

		TClonesArray &prts = *(evt->pip);
		Particle *prt = new(prts[evt->Npip++]) Particle();
		prt->p = pip[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add piminus
	for(unsigned int i=0; i<pim.size(); i++){
		TClonesArray &prts_match = *(evt->pim_match);
		Particle *prt_match = new(prts_match[evt->Npim]) Particle();
		prt_match->p = FindBestMatch(pim[i], pim_match);
		prt_match->x.SetXYZ(0,0,65); // FIXME!!!

		TClonesArray &prts = *(evt->pim);
		Particle *prt = new(prts[evt->Npim++]) Particle();
		prt->p = pim[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}

	// Add proton
	for(unsigned int i=0; i<proton.size(); i++){
		TClonesArray &prts_match = *(evt->proton_match);
		Particle *prt_match = new(prts_match[evt->Nproton]) Particle();
		prt_match->p = FindBestMatch(proton[i], proton_match);
		prt_match->x.SetXYZ(0,0,65); // FIXME!!!

		TClonesArray &prts = *(evt->proton);
		Particle *prt = new(prts[evt->Nproton++]) Particle();
		prt->p = proton[i];
		prt->x.SetXYZ(0,0,65); // FIXME!!!
	}
	

	// Calculate W of reconstructed particles
	for(unsigned int i=0; i<photon.size(); i++)evt->W += photon[i];
	for(unsigned int i=0; i<pip.size(); i++)evt->W += pip[i];
	for(unsigned int i=0; i<pim.size(); i++)evt->W += pim[i];
}

//------------------
// FindBestMatch
//------------------
TLorentzVector DEventProcessor_phys_tree::FindBestMatch(const TLorentzVector &primary, vector<TLorentzVector> &secondaries)
{
	// Loop over secondaries and keep the one with the best figure of merit
	// to return. Initialize return vector with zeros in case not good match
	// is found.
	double max_fom = 0.1;
	TLorentzVector best_match(0.0, 0.0, 0.0);
	for(unsigned int i=0; i<secondaries.size(); i++){
		double fom = GetFOM(primary, secondaries[i]);
		if(fom > max_fom){
			max_fom = fom;
			best_match = secondaries[i];
		}
	}
	
	return best_match;
}

//------------------
// GetFOM
//------------------
double DEventProcessor_phys_tree::GetFOM(const TLorentzVector &a, const TLorentzVector &b) const
{
	// This is a kind of brain-dead algorithm. It wants to use both the
	// momentum direction and magnitude to determine the FOM. For the
	// magnitude, we use the curvature which becomes close to zero for
	// high momentum tracks (good because a 4GeV/c and an 8GeV/c track
	// both look more or less like straight lines).
	//
	// For the direction, we just use the relative angle between the
	// two tracks.
	//
	// For both, we take a reciprocal so that the closer the match,
	// the higher the FOM. We take a product of the 2 FOM components
	// so that both components must have a high value in order for the
	// total FOM to be large.

	double epsilon = 1.0E-6; // prevent reciprocals from resulting in infinity

	double curature_a = 1.0/a.P();
	double curature_b = 1.0/b.P();
	double curvature_diff = fabs(curature_a - curature_b);
	double curvature_fom = 1.0/(curvature_diff + epsilon);

	double theta_rel = fabs(a.Angle(b.Vect()));
	double theta_fom = 1.0/(theta_rel + epsilon);
	
	double fom = curvature_fom*theta_fom;

	return fom;
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

