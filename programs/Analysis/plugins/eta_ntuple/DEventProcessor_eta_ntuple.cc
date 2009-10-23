// $Id: DEventProcessor_eta_ntuple.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_eta_ntuple.cc
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
#include <FCAL/DFCALPhoton.h>
#include <PID/DPhoton.h>
#include <PID/DBeamPhoton.h>

#include "DEventProcessor_eta_ntuple.h"

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
	app->AddProcessor(new DEventProcessor_eta_ntuple());
}
}


//------------------
// init
//------------------
jerror_t DEventProcessor_eta_ntuple::init(void)
{
	// Create tree
	evt = new Event();
	tree = new TTree("t", "#eta Primakoff events");
	tree->Branch("E",&evt);
	
	pthread_mutex_init(&mutex, NULL);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_eta_ntuple::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects
	vector<const DBeamPhoton*> beam_photons;
	vector<const DFCALPhoton*> fcalphotons;
	vector<const DMCThrown*> mcthrowns;

	loop->Get(beam_photons);	// from truth info
	loop->Get(fcalphotons);		// all reconstructed photons in FCAL
	loop->Get(mcthrowns);		// all thrown particles

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Create TLorentzVectors for reconstructed photons
	vector<TLorentzVector> rec_photons;
	vector<TVector3> rec_photons_pos;
	for(unsigned int i=0; i<fcalphotons.size(); i++){
		rec_photons.push_back(fcalphotons[i]->getMom4());
		rec_photons_pos.push_back(fcalphotons[i]->getPosition());
	}

	// Some generators don't supply information on the beam photon. If there
	// are no DBeamPhoton objects, then punt
	if(beam_photons.size()!=1){
		cout<<"Wrong number of DBeamPhoton objects for event "<<eventnumber<<" ("<<beam_photons.size()<<"). Skipping."<<endl;
		return NOERROR;
	}
	TLorentzVector beam_photon = MakeTLorentz(beam_photons[0], 0.0);
	
	// Find thrown eta
	TLorentzVector eta;
	TVector3 vertex;
	bool found_eta=false;
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		if(mcthrowns[i]->pdgtype == 221){
			eta = MakeTLorentz(mcthrowns[i], 0.54745);
			vertex = mcthrowns[i]->position();
			found_eta = true;
			break;
		}
	}
	if(!found_eta){
		cout<<"No thrown eta particle found for event "<<eventnumber<<". Skipping."<<endl;
		return NOERROR;
	}

	pthread_mutex_lock(&mutex);
	
	evt->Clear();
	evt->event = eventnumber;
	evt->beam = beam_photon;
	evt->eta_thrown = eta;
	evt->proton_thrown = target+beam_photon-eta; // assumes coherent!
	evt->vertex = vertex;
	evt->prod_mech = 0;
	evt->decay_mode = 0;
	
	// Loop over reconstructed photons
	for(unsigned int j=0; j<rec_photons.size(); j++){
		evt->AddFCAL(rec_photons[j], rec_photons_pos[j]);
	}
	
	// Loop over all 2-gamma combinations keeping the one closes to the eta mass
	for(unsigned int j=0; j<rec_photons.size(); j++){
		for(unsigned int k=j+1; k<rec_photons.size(); k++){
			if(k==j)continue;
			TLorentzVector my_eta = rec_photons[j] + rec_photons[k];
			if(fabs(evt->eta_best.M()-0.54745)>fabs(my_eta.M()-0.54745))
			evt->eta_best = my_eta;
		}
	}

	// Fill tree
	evt->fcal->Sort(); // sort by cluster energy (uses fcal_t::Compare );
	tree->Fill();

	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// MakeTLorentz
//------------------
TLorentzVector DEventProcessor_eta_ntuple::MakeTLorentz(const DKinematicData *kd, double mass)
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
jerror_t DEventProcessor_eta_ntuple::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_eta_ntuple::fini(void)
{
	// Histograms are automatically written to file by hd_root program (or equivalent)
	// so we don't need to do anything here.

	return NOERROR;
}
