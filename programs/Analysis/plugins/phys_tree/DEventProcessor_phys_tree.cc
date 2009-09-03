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
	event_ptr = &evt;

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
	// Create TRACKING directory
	TDirectory *dir = (TDirectory*)gROOT->FindObject("PHYSICS");
	if(!dir)dir = new TDirectoryFile("PHYSICS","PHYSICS");
	dir->cd();

	tevent = new TTree("event","Event");
	tevent->Branch("E","event",&event_ptr);

	dir->cd("../");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_phys_tree::brun(JEventLoop *loop, int runnumber)
{
	pthread_mutex_lock(&mutex);

	pthread_mutex_unlock(&mutex);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_phys_tree::evnt(JEventLoop *loop, int eventnumber)
{
	// Get reconstructed objects and make TLorentz vectors out of each of them
	vector<const DBeamPhoton*> beam_photons;
	vector<const DPhoton*> photons;
	vector<const DParticle*> particles;
	vector<const DMCThrown*> mcthrowns;

	loop->Get(beam_photons);
	loop->Get(photons);
	loop->Get(particles);
	loop->Get(mcthrowns);

	// Make TLorentzVector for beam photon
	TLorentzVector beam_photon = TLorentzVector(0.0, 0.0, 9.0, 9.0);
	if(beam_photons.size()>0)beam_photon = MakeTLorentz(beam_photons[0], 0.0);
		
	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);

	// Center of mass energy
	double sqrt_s = (beam_photon+target).M();

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
	
	// Total, reconstructed momentum
	TLorentzVector pfinal_recon;
	for(unsigned int i=0; i<rec_piplus.size(); i++)pfinal_recon += rec_piplus[i];
	for(unsigned int i=0; i<rec_piminus.size(); i++)pfinal_recon += rec_piminus[i];
	for(unsigned int i=0; i<rec_protons.size(); i++)pfinal_recon += rec_protons[i];
	
	// Calculate missing momentum
	TLorentzVector missing = beam_photon + target - pfinal_recon;
	
	// If we reconstructed a proton, consider that our proton. Otherwise, assume
	// that's the only thing we missed and use the missing momentum as the proton
	TLorentzVector *proton = rec_protons.size()>0 ? &rec_protons[0]:&missing;
	
	// Mandelstam t
	double t = (*proton - target).Mag2();
	
	// Count the number of thrown charged and neutral particles
	int Nthrown_charged = 0;
	int Nthrown_neutral = 0;
	int Nthrown_charged_fiducial = 0;
	int Nthrown_neutral_fiducial = 0;
	TLorentzVector pthrown_fiducial;
	for(unsigned int k=0; k<mcthrowns.size(); k++){
		TVector3 mom = mcthrowns[k]->momentum();
		if(mcthrowns[k]->charge()==0){
			Nthrown_neutral++;
			if(mom.Mag()<0.100)continue;
			if(mom.Theta()*57.3 < 2.0)continue; // beam hole
			if(mom.Theta()*57.3 > 120)continue; // upstream fiducial limit
			if(fabs(mom.Theta()*57.3-11.0) < 1.0)continue; // BCAL-FCAL gap
			Nthrown_neutral_fiducial++;
			pthrown_fiducial += MakeTLorentz(mcthrowns[k], 0.0);
		}else{
			Nthrown_charged++;
			if(mom.Mag()<0.500)continue;
			if(mom.Theta()*57.3 < 2.0)continue; // beam hole
			if(mom.Theta()*57.3 > 120)continue; // upstream fiducial limit
			Nthrown_charged_fiducial++;
			pthrown_fiducial += MakeTLorentz(mcthrowns[k], mcthrowns[k]->type==14 ? 0.93827:0.13957);
		}
	}
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	evt.event = eventnumber;
	evt.is_fiducial = ((Nthrown_charged+Nthrown_neutral) - (Nthrown_charged_fiducial+Nthrown_neutral_fiducial))<=1;
	evt.Nthrown_charged = Nthrown_charged;
	evt.Nthrown_neutral = Nthrown_neutral;
	evt.Nthrown_charged_fiducial = Nthrown_charged_fiducial;
	evt.Nthrown_neutral_fiducial = Nthrown_neutral_fiducial;
	evt.pthrown_fiducial = pthrown_fiducial;
	evt.pfinal = pfinal_recon;
	evt.pmissing = missing;
	evt.pinitial = beam_photon + target;
	evt.pW = evt.pinitial - *proton;
	evt.W = evt.pW.Mag();
	evt.minus_t = -t;
	evt.sqrt_s = sqrt_s;

	tevent->Fill();

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

