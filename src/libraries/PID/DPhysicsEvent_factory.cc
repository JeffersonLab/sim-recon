// $Id$
//
//    File: DPhysicsEvent_factory.cc
// Created: Wed Aug  4 10:37:55 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.4.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <PID/DChargedTrack.h>
#include <PID/DPhoton.h>
#include <PID/DVertex.h>

#include "DPhysicsEvent_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPhysicsEvent_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPhysicsEvent_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	// This is not a strict limit. Rather, it is the size the pool will be reduced
	// to if it has grown larger on the previous event and the current event does
	// not require this many. This prevents memory-leak-like behavior when running
	// many threads where each thread keeps allocating bigger and bigger pools as
	// it comes across slightly busier and busier events.
	MAX_PARTINFOS = 10;

	Nbinst = 600;
	tmin = -100.0;
	tmax = 500.0;
	Nbinsz = 50;
	zmin = 0.0;
	zmax = 100.0;
	
	MAKE_ROOT_HISTOS = false;
	gPARMS->SetDefaultParameter("PHYSICS:MAKE_ROOT_HISTOS", MAKE_ROOT_HISTOS);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPhysicsEvent_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Create a 2-D histogram of time vs. z-position for each particle
	// at the vertex position. All histos will be added and any maximum
	// found in the sum will be used to identify particles belonging
	// to the group. Remaining particles will have their histos added
	// and the process repeated until all groups are found.

	// Get reconstructed charged particles and photons
	vector<const DChargedTrack*> chargedtracks;
	vector<const DPhoton*> all_photons;
	loop->Get(chargedtracks);
	loop->Get(all_photons);

	// Keep only most probable hypothesis for charged tracks
	vector<const DTrackTimeBased*> tracktimebaseds;
	for(unsigned int i=0; i<chargedtracks.size(); i++){
		const vector<const DTrackTimeBased*> &hypotheses = chargedtracks[i]->hypotheses;
		if(hypotheses.size()>0)tracktimebaseds.push_back(hypotheses[0]);
	}

	// Keep only photons not matched to charged tracks
	vector<const DPhoton*> photons;
	for(unsigned int i=0; i<all_photons.size(); i++){
		if(all_photons[i]->getTag()!=DPhoton::kCharge)photons.push_back(all_photons[i]);
	}

	// To minimize memory usage and time in allocation, we maintain a
	// pool of partInfo_t objects. Make sure the pool is large enough to hold
	// all of the particles we have this event. 
	unsigned int Nparticles_total = tracktimebaseds.size() + photons.size();
	for(unsigned int i=partInfos_pool.size(); i<Nparticles_total; i++){
		partInfo_t *pi = new partInfo_t();

		pi->SetLimits(tmin, tmax, zmin, zmax, Nbinst, Nbinsz);
		partInfos_pool.push_back(pi);
	}

	// Periodically delete some partInfo_t objects if the pool gets too large.
	// This prevents memory-leakage-like behavor.
	if((Nparticles_total < MAX_PARTINFOS) && (partInfos_pool.size()>MAX_PARTINFOS)){
		for(unsigned int i=MAX_PARTINFOS; i<partInfos_pool.size(); i++)delete partInfos_pool[i];
		partInfos_pool.resize(MAX_PARTINFOS);
	}
	
	// Vector to hold list of partInfo_t objects for all particles (photons
	// and charged tracks)
	vector<partInfo_t*> parts;

	// Assign and fill partInfo_t objects for each charged track
	for(unsigned int i=0; i<tracktimebaseds.size(); i++){
		partInfo_t *pi = partInfos_pool[parts.size()];
		FillPartInfoChargedTrack(pi, tracktimebaseds[i]);
		parts.push_back(pi);
	}

	// Assign and fill partInfo_t objects for each photon
	for(unsigned int i=0; i<photons.size(); i++){
		partInfo_t *pi = partInfos_pool[parts.size()];
		FillPartInfoPhoton(pi, photons[i]);
		parts.push_back(pi);
	}
	
	// At this point we have combined all photons and charged tracks into
	// a single list of objects (parts). Each has a histogram of t vs.z
	// values filled using approriate uncertainties (no covariance). We
	// can now use this list to identify resonances in the t/z plane which
	// indicate a vertex location. Particles within 3sigma in both t and
	// z of the resonance will be grouped together as belonging to the 
	// same physics event. We loop until all particles have been assigned
	// to a group, even if that means assigning particles to their own
	// "group of one".

	// Loop until all particles have been assigned to a group.
	vector< vector<partInfo_t *> > groups;
	do{
		// Make a list of all particles that have not been assigned
		// to a group
		vector<const DHoughFind*> unassigned;
		for(unsigned int i=0; i<parts.size(); i++){
			if(!parts[i]->is_in_group)unassigned.push_back(parts[i]);
		}
		
		// Find the maximum t,z coordinate by adding all unassigned
		// particle's histos together
		DVector2 maxloc = DHoughFind::GetMaxBinLocation(unassigned);
		
		if(debug_level>0)_DBG_<<"Location of maximum: t="<<maxloc.X()<<"  z="<<maxloc.Y()<<endl;		

		// Loop over all unassigned particles, assigning any within
		// 3 sigma in both t and z to the new group. We loop over
		// the parts vector just because it saves a dynamic_cast
		// if we were to use the unassigned vector.
		vector<partInfo_t *> new_group;
		for(unsigned int i=0; i<parts.size(); i++){
			partInfo_t *pi = parts[i];
			if(pi->is_in_group)continue;
			
			double delta_t = fabs(maxloc.X() - pi->t);
			if(delta_t/pi->sigmat > 3.0)continue;
			
			double delta_z = fabs(maxloc.Y() - pi->z);
			if(delta_z/pi->sigmaz > 3.0)continue;

			// Assign this particle to the group
			new_group.push_back(pi);
		}
		
		// At this point it's possible (but hopefully unlikely) that the
		// maximum in the t,z sum histo was generated at an in-between place
		// with no single particle nearby. In that case, the new_group is
		// empty, even though there are unassigned particles. The best we
		// can do here is to assign one particle to the new_group and hope
		// that the next iteration groups the remaining ones appropriately.
		// To try and minimize the chances of placing a particle from the
		// L1 trigger event in its own group, we choose the particle with a time
		// the furthest away from t=0.
		if(new_group.size()==0){
			partInfo_t *pi_with_max_t = NULL;
			double delta_t_max=0.0;
			for(unsigned int i=0; i<parts.size(); i++){
				partInfo_t *pi = parts[i];
				if(pi->is_in_group)continue;
			
				double delta_t = fabs(maxloc.X() - pi->t);
				if(delta_t>delta_t_max || pi_with_max_t==NULL){
					delta_t_max = delta_t;
					pi_with_max_t = pi;
				}
			}
			
			if(pi_with_max_t==NULL){
				_DBG_<<"pi_with_max_t==NULL. This should never happen! Complain to davidl@jlab.org"<<endl;
				break;
			}
			
			new_group.push_back(pi_with_max_t);
		}
		
		// Set the is_in_group flags for all of the members of the new group
		for(unsigned int i=0; i<new_group.size(); i++)new_group[i]->is_in_group = true;
		
		// Add the new group to the list of groups
		groups.push_back(new_group);
		
	}while(!AllInGroups(parts));

	// OK, we've now grouped the particles together into groups. Create a new
	// DPhysicsEvent for each group
	for(unsigned int i=0; i<groups.size(); i++){
		vector<partInfo_t *> &group = groups[i];
		
		DPhysicsEvent *pe = new DPhysicsEvent;
		
		// Loop over particles in the group, pushing them onto the
		// appropriate vectors in the DPhysicsEvent
		for(unsigned int j=0; j<group.size(); j++){
			partInfo_t *pi = group[j];
			
			if(pi->photon){
				// photon
				pe->photon.push_back(pi->photon);
			}else if(pi->track){
				// Determine type of charged track
				double mass = pi->track->mass();
				double q = pi->track->charge();
				if(fabs(mass-0.13957)<0.001){
					if(q>0.0){
						// pi+
						pe->pip.push_back(pi->track);
					} else {
						// pi-
						pe->pim.push_back(pi->track);
					}
				} else if(fabs(mass-0.93827)<0.001){
					// proton
					pe->proton.push_back(pi->track);
				} else if(fabs(mass-0.49368)<0.001){
					if(q>0.0){
						// K+
						pe->Kp.push_back(pi->track);
					} else {
						// K-
						pe->Km.push_back(pi->track);
					}
				} else {
					if(q>0.0){
						// other +
						pe->otherp.push_back(pi->track);
					} else {
						// other -
						pe->otherm.push_back(pi->track);
					}
				}
			}else{
				_DBG_<<"Both photon and track points in partInfo_t are NULL. This should never happen. Compalin to davidl@jlab.org"<<endl;
			}
		}
		
		// "Publish" the DPhysicsEvent object 
		_data.push_back(pe);
	}
	
	// Optionally record histo info for debugging
	if(MAKE_ROOT_HISTOS)MakeRootHists(eventnumber, groups);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPhysicsEvent_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPhysicsEvent_factory::fini(void)
{
	return NOERROR;
}

//------------------
// AllInGroups
//------------------
bool DPhysicsEvent_factory::AllInGroups(vector<partInfo_t*> &parts)
{
	for(unsigned int i=0; i<parts.size(); i++)if(!parts[i]->is_in_group)return false;
	
	return true;
}

//------------------
// FillPartInfoChargedTrack
//------------------
void DPhysicsEvent_factory::FillPartInfoChargedTrack(DPhysicsEvent_factory::partInfo_t *pi, const DTrackTimeBased *trk)
{
	pi->Reset();
	pi->track = trk;

	pi->t = trk->t0();
	pi->sigmat = trk->t0_err();
	pi->z = trk->z();
	pi->sigmaz = 0.8/sin(trk->momentum().Theta()); // in cm.  For now, use 3mm wide angle track resolution scaled by sin(theta)

	pi->Fill(pi->t, pi->sigmat, pi->z, pi->sigmaz);
}

//------------------
// FillPartInfoPhoton
//------------------
void DPhysicsEvent_factory::FillPartInfoPhoton(DPhysicsEvent_factory::partInfo_t *pi, const DPhoton *photon)
{
	pi->Reset();
	pi->photon = photon;

	pi->t = photon->t0();
	pi->sigmat = photon->t0_err();
	pi->z = photon->z();
	pi->sigmaz = 30.0/sqrt(12); // in cm. Use length of target for z-resolution of photons

	pi->Fill(pi->t, pi->sigmat, pi->z, pi->sigmaz);
}

//------------------
// MakeRootHists
//------------------
void DPhysicsEvent_factory::MakeRootHists(int event, vector< vector<partInfo_t *> > &groups)
{
	/// This is meant for debugging ONLY. It will take the DHoughFind objects
	/// contained in the given "groups" container and convert them into
	/// ROOT histograms. This should be utilized using something like hd_root
	/// so that the histograms are saved to a file.
	///
	/// A TDirectory will be made to hold the histograms for each event.
	/// Within that, a TDirectory will be created for each group and the
	/// group member's histos will be created there.

	// Save current directory so we can cd back to it before returing
	TDirectory *saveDir = gDirectory;

	// Create a TDirectory to hold the event
	char dirname[256];
	sprintf(dirname, "event%03d", event);
	TDirectory *eventdir = saveDir->mkdir(dirname);

	// Container to keep pointer to all hists so we can make sum histo later
	vector<TH2D*> root_hists;

	// Loop over groups
	for(unsigned int i=0; i<groups.size(); i++){
		
		// Make a TDirectory for this group
		sprintf(dirname, "group%02d", i);
		TDirectory *groupdir = eventdir->mkdir(dirname);
		groupdir->cd();
		
		// Loop over members of group
		vector<partInfo_t *> &group = groups[i];
		for(unsigned int j=0; j<group.size(); j++){
			char hname[256];
			sprintf(hname, "part%02d", j);
			TH2D *h = group[j]->MakeIntoRootHist(hname);
			h->SetDrawOption("surf4");
			h->SetXTitle("time from L1 trigger (ns)");
			h->SetYTitle("Z position of vertex (cm)");
			h->SetStats(0);
			root_hists.push_back(h);
		}
	}

	// Create sum histo for all particles this event
	if(root_hists.size() > 0){
		eventdir->cd();
		TH2D *sum = (TH2D*)root_hists[0]->Clone("sum");
		for(unsigned int i=1; i<root_hists.size(); i++){
			sum->Add(root_hists[i]);
		}
	}

	// Write eventdir and its contents to the file
	//eventdir->Write();

	// cd back into TDirectory we were in upon entering this method
	saveDir->cd();
}


