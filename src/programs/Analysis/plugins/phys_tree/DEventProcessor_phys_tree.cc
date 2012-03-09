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
#include <PID/DBeamPhoton.h>
#include <PID/DParticleSet.h>
#include <PID/DPhysicsEvent.h>
#include <TRACKING/DMCThrown.h>

bool static CompareLorentzEnergy(const Particle &a, const Particle &b){
  return a.p.E()<b.p.E();
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
	vector<const DPhysicsEvent*> physicsevents;

	loop->Get(beam_photons);
	loop->Get(mcthrowns);
	loop->Get(physicsevents);

	TVector3 VertexRec, VertexGen;
	VertexGen = TVector3(mcthrowns[0]->position().X(),mcthrowns[0]->position().Y(),mcthrowns[0]->position().Z());
	// Make Particle object for beam photon
	TLorentzVector beam_photon(0.0, 0.0, 9.0, 9.0);
	if(beam_photons.size()>0){
		const DLorentzVector &lv = beam_photons[0]->lorentzMomentum();
		beam_photon.SetPxPyPzE(lv.Px(), lv.Py(), lv.Pz(), lv.E());
	}

	// Target is proton at rest in lab frame
	TLorentzVector target(0.0, 0.0, 0.0, 0.93827);
	
	// Find the DPhysicsEvent with the most particles and only use that one.
	// This is not a long term solution, but is motivated by the fact that
	// we have only one set of DMCThrown particles and one DBeamPhoton
	const DPhysicsEvent *physicsevent = NULL;
	int max_parts=0;
	for(unsigned int i=0; i<physicsevents.size(); i++){
		const DParticleSet *ps = physicsevents[i]->particle_sets[0];
		int Nparts = ps->pip.size() + ps->pim.size()
		        + ps->photon.size() + ps->proton.size()
		        + ps->Kp.size()     + ps->Km.size()
		        + ps->otherp.size() + ps->othern.size()
		        + ps->otherz.size() + ps->neutron.size();
		if(Nparts>max_parts || physicsevents[i]==NULL){
			physicsevent = physicsevents[i];
			max_parts = Nparts;
		}
		VertexRec.SetX(ps->vertex->dSpacetimeVertex.X());
		VertexRec.SetY(ps->vertex->dSpacetimeVertex.Y());
		VertexRec.SetZ(ps->vertex->dSpacetimeVertex.Z());
	}

	// Create Particle objects for each of the common particle types
	particle_set rec;

	if (physicsevent!=NULL){
	  const DParticleSet *particle_set=physicsevent->particle_sets[0];
	  for(unsigned int j=0; j<particle_set->photon.size(); j++){
		// photon
	    rec.photons.push_back(MakeParticle(particle_set->photon[j]->dNeutralParticleHypotheses[0], 0.0));
	  }
	  for(unsigned int j=0; j<particle_set->neutron.size(); j++){
		// neutron
	    rec.neutrons.push_back(MakeParticle(particle_set->neutron[j]->dNeutralParticleHypotheses[0], 0.939565));
	  }
	  for(unsigned int j=0; j<particle_set->pip.size(); j++){
	    // pi+
	    rec.piplus.push_back(MakeParticle(particle_set->pip[j]->dChargedTrackHypotheses[0], 0.13957));
	  }
	  for(unsigned int j=0; j<particle_set->pim.size(); j++){
	    // pi-
	    rec.piminus.push_back(MakeParticle(particle_set->pim[j]->dChargedTrackHypotheses[0], 0.13957));
	  }
	  for(unsigned int j=0; j<particle_set->proton.size(); j++){
	    // proton
	    rec.protons.push_back(MakeParticle(particle_set->proton[j]->dChargedTrackHypotheses[0], 0.93827));
	  }
	  for(unsigned int j=0; j<particle_set->Kp.size(); j++){
	    // K+
	    rec.Kplus.push_back(MakeParticle(particle_set->Kp[j]->dChargedTrackHypotheses[0], 0.493677));
	  }
	  for(unsigned int j=0; j<particle_set->Km.size(); j++){
	    // K-
	    rec.Kminus.push_back(MakeParticle(particle_set->Km[j]->dChargedTrackHypotheses[0], 0.493677));
	  }  
	}

	// Create Particle objects for thrown particles
	bool all_mesons_fiducial = true;
	bool all_photons_fiducial = true;
	bool all_neutrons_fiducial = true;
	bool all_protons_fiducial = true;
	particle_set thr;
	for(unsigned int k=0; k<mcthrowns.size(); k++){
	  switch(mcthrowns[k]->type){
	  case  1: thr.photons.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.0));
	    all_photons_fiducial &= IsFiducial(mcthrowns[k]);
	    break;
	  case  8: thr.piplus.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.13957));
	      all_mesons_fiducial &= IsFiducial(mcthrowns[k]);
	      break;
	  case  9: thr.piminus.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.13957));
	    all_mesons_fiducial &= IsFiducial(mcthrowns[k]);
	    break;
	  case 11: thr.Kplus.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.493677));
	    all_mesons_fiducial &= IsFiducial(mcthrowns[k]);
	    break;
	  case 12: thr.Kminus.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.493677));
	    all_mesons_fiducial &= IsFiducial(mcthrowns[k]);
	    break;
	  case 13: thr.neutrons.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.939565));
	    all_neutrons_fiducial &= IsFiducial(mcthrowns[k]);
	    break;
	  case 14: thr.protons.push_back(MakeParticle((DKinematicData*)mcthrowns[k], 0.93827));
	    all_protons_fiducial &= IsFiducial(mcthrowns[k]);
	      break;
	  }
	}
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Fill in Event objects for both thrown and reconstructed
	evt_recon->Clear();
	evt_thrown->Clear();
	FillEvent(evt_recon, rec, thr);
	FillEvent(evt_thrown, thr, rec);

	// Copy fiducial cuts (based only on thrown values) to both trees
	bool all_fiducial = all_mesons_fiducial && all_photons_fiducial && all_protons_fiducial && all_neutrons_fiducial;
	evt_recon->all_fiducial = all_fiducial;
	evt_recon->all_mesons_fiducial = all_mesons_fiducial;
	evt_recon->all_photons_fiducial = all_photons_fiducial;
	evt_recon->all_neutrons_fiducial = all_neutrons_fiducial;
	evt_recon->all_protons_fiducial = all_protons_fiducial;
	evt_recon->vertex = VertexRec;

	evt_thrown->all_fiducial = all_fiducial;
	evt_thrown->all_mesons_fiducial = all_mesons_fiducial;
	evt_thrown->all_photons_fiducial = all_photons_fiducial;
	evt_thrown->all_neutrons_fiducial = all_neutrons_fiducial;
	evt_thrown->all_protons_fiducial = all_protons_fiducial;
	evt_thrown->vertex = VertexGen;

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
// MakeParticle
//------------------
Particle DEventProcessor_phys_tree::MakeParticle(const DKinematicData *kd, double mass)
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
	double x = kd->position().X();
	double y = kd->position().Y();
	double z = kd->position().Z();
	
	Particle part;
	part.p.SetPxPyPzE(px,py,pz,E);
	part.x.SetXYZ(x, y, z);
	part.is_fiducial = IsFiducial(kd);
	part.chisq = -1.0;
	part.Ndof = 0;
	part.FOM_pid = -1.0;
	
	return part;
}

//------------------
// MakeParticle
//------------------
Particle DEventProcessor_phys_tree::MakeParticle(const DChargedTrackHypothesis *locChargedTrackHypothesis, double mass)
{
	// Most values get set using DKinematicData part
	Particle part = MakeParticle((DKinematicData*)locChargedTrackHypothesis, mass);

	// Things specific to DChargedTrackHypothesis
	part.chisq = locChargedTrackHypothesis->dChiSq;
	part.Ndof = locChargedTrackHypothesis->dNDF;
	part.FOM_pid = locChargedTrackHypothesis->dFOM;

	return part;
}

//------------------
// MakeParticle
//------------------
Particle DEventProcessor_phys_tree::MakeParticle(const DNeutralParticleHypothesis *locNeutralParticleHypothesis, double mass)
{
	// Most values get set using DKinematicData part
	Particle part = MakeParticle((DKinematicData*)locNeutralParticleHypothesis, mass);

	// Things specific to DNeutralParticleHypothesis
	part.chisq = locNeutralParticleHypothesis->dChiSq;
	part.Ndof = locNeutralParticleHypothesis->dNDF;
	part.FOM_pid = locNeutralParticleHypothesis->dFOM;

	return part;
}

//------------------
// FillEvent
//------------------
void DEventProcessor_phys_tree::FillEvent(Event *evt, particle_set &pset, particle_set &pset_match)
{

	vector<Particle> &photon = pset.photons;
	vector<Particle> &neutron = pset.neutrons;
	vector<Particle> &pip = pset.piplus;
	vector<Particle> &pim = pset.piminus;
	vector<Particle> &Kp = pset.Kplus;
	vector<Particle> &Km = pset.Kminus;
	vector<Particle> &proton = pset.protons;

	vector<Particle> &photon_match = pset_match.photons;
	vector<Particle> &neutron_match = pset_match.neutrons;
	vector<Particle> &pip_match = pset_match.piplus;
	vector<Particle> &pim_match = pset_match.piminus;
	vector<Particle> &Kp_match = pset_match.Kplus;
	vector<Particle> &Km_match = pset_match.Kminus;
	vector<Particle> &proton_match = pset_match.protons;

	// Sort particle arrays by energy
	sort(photon.begin(), photon.end(), CompareLorentzEnergy);
	sort(neutron.begin(), neutron.end(), CompareLorentzEnergy);
	sort(pip.begin(), pip.end(), CompareLorentzEnergy);
	sort(pim.begin(), pim.end(), CompareLorentzEnergy);
	sort(Kp.begin(), Kp.end(), CompareLorentzEnergy);
	sort(Km.begin(), Km.end(), CompareLorentzEnergy);
	sort(proton.begin(), proton.end(), CompareLorentzEnergy);

	// Add photons
	for(unsigned int i=0; i<photon.size(); i++){
		TClonesArray &prts_match = *(evt->photon_match);
		Particle *prt_match = new(prts_match[evt->Nphoton]) Particle();
		*prt_match = FindBestMatch(photon[i], photon_match);

		TClonesArray &prts = *(evt->photon);
		Particle *prt = new(prts[evt->Nphoton++]) Particle();
		*prt = photon[i];
	}

	// Add neutrons
	for(unsigned int i=0; i<neutron.size(); i++){
		TClonesArray &prts_match = *(evt->neutron_match);
		Particle *prt_match = new(prts_match[evt->Nneutron]) Particle();
		*prt_match = FindBestMatch(neutron[i], neutron_match);

		TClonesArray &prts = *(evt->neutron);
		Particle *prt = new(prts[evt->Nneutron++]) Particle();
		*prt = neutron[i];
	}

	// Add piplus
	for(unsigned int i=0; i<pip.size(); i++){
		TClonesArray &prts_match = *(evt->pip_match);
		Particle *prt_match = new(prts_match[evt->Npip]) Particle();
		*prt_match = FindBestMatch(pip[i], pip_match);

		TClonesArray &prts = *(evt->pip);
		Particle *prt = new(prts[evt->Npip++]) Particle();
		*prt = pip[i];
	}

	// Add piminus
	for(unsigned int i=0; i<pim.size(); i++){
		TClonesArray &prts_match = *(evt->pim_match);
		Particle *prt_match = new(prts_match[evt->Npim]) Particle();
		*prt_match = FindBestMatch(pim[i], pim_match);

		TClonesArray &prts = *(evt->pim);
		Particle *prt = new(prts[evt->Npim++]) Particle();
		*prt = pim[i];
	}

	// Add Kplus
	for(unsigned int i=0; i<Kp.size(); i++){
		TClonesArray &prts_match = *(evt->Kp_match);
		Particle *prt_match = new(prts_match[evt->NKp]) Particle();
		*prt_match = FindBestMatch(Kp[i], Kp_match);

		TClonesArray &prts = *(evt->Kp);
		Particle *prt = new(prts[evt->NKp++]) Particle();
		*prt = Kp[i];
	}

	// Add Kminus
	for(unsigned int i=0; i<Km.size(); i++){
		TClonesArray &prts_match = *(evt->Km_match);
		Particle *prt_match = new(prts_match[evt->NKm]) Particle();
		*prt_match = FindBestMatch(Km[i], Km_match);

		TClonesArray &prts = *(evt->Km);
		Particle *prt = new(prts[evt->NKm++]) Particle();
		*prt = Km[i];
	}

	// Add proton
	for(unsigned int i=0; i<proton.size(); i++){
		TClonesArray &prts_match = *(evt->proton_match);
		Particle *prt_match = new(prts_match[evt->Nproton]) Particle();
		*prt_match = FindBestMatch(proton[i], proton_match);

		TClonesArray &prts = *(evt->proton);
		Particle *prt = new(prts[evt->Nproton++]) Particle();
		*prt = proton[i];
	}
	
	// Calculate W of reconstructed particles
	for(unsigned int i=0; i<photon.size(); i++)evt->W += photon[i].p;
	for(unsigned int i=0; i<neutron.size(); i++)evt->W += neutron[i].p;
	for(unsigned int i=0; i<pip.size(); i++)evt->W += pip[i].p;
	for(unsigned int i=0; i<pim.size(); i++)evt->W += pim[i].p;
	for(unsigned int i=0; i<Kp.size(); i++)evt->W += Kp[i].p;
	for(unsigned int i=0; i<Km.size(); i++)evt->W += Km[i].p;
}

//------------------
// FindBestMatch
//------------------
Particle DEventProcessor_phys_tree::FindBestMatch(const Particle &primary, vector<Particle> &secondaries)
{
	// Loop over secondaries and keep the one with the best figure of merit
	// to return. Initialize return vector with zeros in case not good match
	// is found.
	double max_fom = 0.1;
	Particle best_match;
	best_match.p.SetXYZT(0.0, 0.0, 0.0, 0.0);
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
double DEventProcessor_phys_tree::GetFOM(const Particle &a, const Particle &b) const
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

	double curature_a = 1.0/a.p.P();
	double curature_b = 1.0/b.p.P();
	double curvature_diff = fabs(curature_a - curature_b);
	double curvature_fom = 1.0/(curvature_diff + epsilon);

	double theta_rel = fabs(a.p.Angle(b.p.Vect()));
	double theta_fom = 1.0/(theta_rel + epsilon);
	
	double fom = curvature_fom*theta_fom;

	return fom;
}

//------------------
// IsFiducial
//------------------
bool DEventProcessor_phys_tree::IsFiducial(const DKinematicData *kd)
{
	double theta_degrees = kd->momentum().Theta()*TMath::RadToDeg();
	double p = kd->momentum().Mag();
	
	if(kd->charge()==0.0){
		// photon
		if(theta_degrees<2.0 || theta_degrees>110.0)return false;
		if(p<0.100)return false;
	}else{
		// charged particle
		if(theta_degrees<1.0 || theta_degrees>120.0)return false;
		if(p<0.400)return false;
	}
	
	return true;
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

