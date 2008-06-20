// $Id$
//
//    File: DEventProcessor_track_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_track_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrack.h>
#include <FDC/DFDCGeometry.h>
#include <DVector2.h>
#include <particleType.h>


// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_track_hists());
}
} // "C"


//------------------
// DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::DEventProcessor_track_hists()
{	
	trk_ptr = &trk;
	cdchit_ptr = &cdchit;
	ref=NULL;
	
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&rt_mutex, NULL);
}

//------------------
// ~DEventProcessor_track_hists
//------------------
DEventProcessor_track_hists::~DEventProcessor_track_hists()
{
	if(ref)delete ref;
}

//------------------
// init
//------------------
jerror_t DEventProcessor_track_hists::init(void)
{
	// Create Trees
	trkeff = new TTree("trkeff","Tracking Efficiency");
	trkeff->Branch("F","track",&trk_ptr);

	cdchits = new TTree("cdchits","CDC hits");
	cdchits->Branch("C","dchit",&cdchit_ptr);
	
	MAX_HIT_DIST_CDC = 1.0; // cm

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_track_hists::brun(JEventLoop *loop, int runnumber)
{
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	bfield = dapp->GetBfield(); // temporary until new geometry scheme is worked out
	ref = new DReferenceTrajectory(bfield);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_track_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_track_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_track_hists::evnt(JEventLoop *loop, int eventnumber)
{
	CDChitv cdctrackhits;
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DMCThrown*> mcthrowns;
	
	loop->Get(cdctrackhits);
	loop->Get(trackcandidates, "CDC");
	loop->Get(mcthrowns);
	
	// Lock mutex
	pthread_mutex_lock(&mutex);

	// Fill maps associating CDC and FDC hits with truth points
	FindCDCTrackNumbers(loop);
		
	// Get hit list for all candidates
	vector<CDChitv> cdc_candidate_hits;
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		CDChitv cdc_outhits;
		GetCDCHits(trackcandidates[i], cdctrackhits, cdc_outhits);
		cdc_candidate_hits.push_back(cdc_outhits);
	}

	// Get hit list for all throwns
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// if this isn't a charged track, then skip it
		if(fabs(mcthrowns[i]->charge())==0.0)continue;

		CDChitv cdc_thrownhits;
		GetCDCHitsFromTruth(i+1, cdc_thrownhits);

		trk.status_can = 0;
		if(cdc_thrownhits.size()<5){
			trk.status_can--;
		}
		
		trk.pthrown = mcthrown->momentum();
		trk.z_thrown = mcthrown->position().Z();
		trk.ncdc_hits_thrown = cdc_thrownhits.size();
		trk.ncdc_hits = cdctrackhits.size();
		
		//========= CANDIDATE ========
		// Look for a candidate track that matches this thrown one
		CDChitv cdc_matched_hits;
		unsigned int icdc_can = FindMatch(cdc_thrownhits, cdc_candidate_hits, cdc_matched_hits);
		
		// Initialize values for when no matching candidate is found
		trk.pcan.SetXYZ(0.0, 0.0, 0.0);
		trk.z_can = -1000.0;
		trk.ncdc_hits_can = 0;
		trk.ncdc_hits_thrown_and_can = 0;
		trk.cdc_chisq_can = 1.0E6;
		
		// CDC
		const DTrackCandidate *cdc_can = NULL;
		if(icdc_can>=0 && icdc_can<trackcandidates.size()){
			cdc_can = trackcandidates[icdc_can];
			
			trk.ncdc_hits_can = cdc_candidate_hits[icdc_can].size();
			trk.ncdc_hits_thrown_and_can = cdc_matched_hits.size();
			trk.cdc_chisq_can = 0.0;
		}
		
		// Figure out which candidate (if any) I should match this with
		const DTrackCandidate* can = cdc_can;
		
		if(can!=NULL){
			trk.pcan = can->momentum();
			trk.z_can = can->position().Z();
		}else{
			trk.pcan.SetXYZ(0,0,0);
			trk.z_can = -1000.0;
		}

		
		// Fill tree
		trkeff->Fill();
	}
	
	// CDC hits tree
	map<const DCDCTrackHit*, const DMCTrackHit*>::iterator iter;
	for(iter=cdclink.begin(); iter!=cdclink.end(); iter++){
		const DMCTrackHit *truth = iter->second;
		const DCDCTrackHit *hit = iter->first;
		const DCDCWire *wire = hit->wire;

		// Fill in values we can with only hit and truth point
		TVector3 pos_truth(truth->r*cos(truth->phi), truth->r*sin(truth->phi), truth->z);
		cdchit.wire = wire->straw;
		cdchit.layer = wire->ring;
		cdchit.t = hit->tdrift;
		cdchit.pos_truth = pos_truth;
		cdchit.ptype = truth->ptype;

		// Some quantities require a fit track to calculate (e.g. residuals)
		// Here we see if there is a fit track associated with the thrown
		// track this hit came from.
		int trackno = truth->track;
		const DMCThrown *mcthrown=NULL;
		if(trackno>=1 && trackno<=(int)mcthrowns.size()){
			mcthrown=mcthrowns[trackno-1];

			// We just need the thrown track for the beta
			double mass = ParticleMass((Particle_t)mcthrown->type);
			cdchit.beta = sqrt(1.0/(pow(mass/mcthrown->momentum().Mag(), 2.0) + 1.0));
		}else{
			cdchit.beta = -1000.0;
		}
		
		// Now, fill in the rest of the DC hit info based on whether we have
		// a fit track or not

		// No fit track
		cdchit.tof = -1000.0;
		cdchit.doca = -1000.0;
		cdchit.resi = -1000.0;
		cdchit.track_wire_angle = -1000.0;
		cdchit.chisq = -1000.0;
		cdchit.pos_wire.SetXYZ(-1000,-1000,-1000);
		cdchit.pos_doca.SetXYZ(-1000,-1000,-1000);
		
		cdchits->Fill();
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}


//------------------
// GetCDCHits
//------------------
void DEventProcessor_track_hists::GetCDCHits(const DKinematicData *p, CDChitv &inhits, CDChitv &outhits)
{
	// In case we run with multiple threads
	pthread_mutex_lock(&rt_mutex);
	
	// Re-swim the reference trajectory using the parameters in p
	ref->Swim(p->position(), p->momentum(), p->charge());
	
	// Loop over hits in the "inhits" vector and find the DOCA of the R.T.
	// for each.
	outhits.clear();
	if(ref->Nswim_steps == 0){
		pthread_mutex_unlock(&rt_mutex);
		return;
	}
	for(unsigned int i=0; i<inhits.size(); i++){
		double doca = ref->DistToRT(inhits[i]->wire);
		if(doca < MAX_HIT_DIST_CDC)outhits.push_back(inhits[i]);
	}
	
	pthread_mutex_unlock(&rt_mutex);
}

//------------------
// GetCDCHitsFromTruth
//------------------
void DEventProcessor_track_hists::GetCDCHitsFromTruth(int trackno, CDChitv &outhits)
{
	// Loop over all entries in the cdclink map and find the ones
	// corresponding to the given track number
	outhits.clear();
	map<const DCDCTrackHit*, const DMCTrackHit*>::iterator iter;
	for(iter=cdclink.begin(); iter!=cdclink.end(); iter++){
		if((iter->second)->track==trackno)outhits.push_back(iter->first);
	}
}

//------------------
// FindMatch
//------------------
unsigned int DEventProcessor_track_hists::FindMatch(CDChitv &thrownhits, vector<CDChitv> &candidate_hits, CDChitv &matched_hits)
{
	// Loop through all of the candidate hits and look for the one that best
	// matches this thrown's hits. The algorithm simply looks for the candidate
	// that has the most hits in common with the thrown.
	//
	// If the best match shares less than half its hits with the thrown,
	// then it is considered not to match and -1 is returned.

	unsigned int ibest=(unsigned int)-1;
	CDChitv my_matched;
	matched_hits.clear();
	for(unsigned int i=0; i<candidate_hits.size(); i++){
		CDChitv &canhits = candidate_hits[i];
		
		// Loop over thrown hits
		my_matched.clear();
		for(unsigned int j=0; j<thrownhits.size(); j++){
			// Loop over candidate hits
			for(unsigned int k=0; k<canhits.size(); k++){
				if(canhits[k] == thrownhits[j])my_matched.push_back(thrownhits[j]);
			}
		}
		
		// Check if this candidate is a better match
		if(my_matched.size() > matched_hits.size()){
			matched_hits = my_matched;
			ibest = i;
		}
	}
	
	// Is the best good enough?
	if(matched_hits.size() >= (thrownhits.size()/2))return ibest;
	
	return (unsigned int)-1;
}

//------------------
// FindCDCTrackNumbers
//------------------
void DEventProcessor_track_hists::FindCDCTrackNumbers(JEventLoop *loop)
{
	vector<const DCDCTrackHit*> cdchits;
	vector<const DMCTrackHit*> mchits;
	loop->Get(cdchits);
	loop->Get(mchits);

	cdclink.clear();

	// Loop over all CDC wire hits
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCTrackHit *cdchit = cdchits[i];
		const DCDCWire *wire = cdchit->wire;
		if(!wire)continue;
		
		// Loop over CDC truth points
		for(unsigned int j=0; j<mchits.size(); j++){
			const DMCTrackHit *mchit = mchits[j];
			if(mchit->system != SYS_CDC)continue;
			
			// Find the distance between this truth hit and this wire
			DVector3 pos_truth(mchit->r*cos(mchit->phi), mchit->r*sin(mchit->phi), mchit->z);
			DVector3 A = wire->udir.Cross(pos_truth - wire->origin);
			double dist = A.Mag();
			
			// The value of cdchit->dist was calculated using a time
			// that did not have TOF removed. We must estimate the TOF
			// here so so we can correct for that in order to make a
			// more accurate comparison between the truth point
			// and the hit. We estimate TOF by finding the distance
			// that the truth point is from the center of the target
			// and assume beta=1.
			DVector3 target(0,0,65.0);
			DVector3 delta_truth = pos_truth - target;
			double tof = delta_truth.Mag()/3.0E10/1.0E-9;
			double corrected_dist = cdchit->dist*(cdchit->tdrift-tof)/cdchit->tdrift;
//_DBG_<<"x="<<pos_truth.X()<<" y="<<pos_truth.Y()<<" z="<<pos_truth.Z()<<" tof="<<tof<<" tdrift="<<cdchit->tdrift-tof<<" t="<<cdchit->tdrift<<endl;
			// The best we can do is check that the truth hit is inside
			// the straw.
			bool match = fabs(corrected_dist - dist)<0.8;
//_DBG_<<"fabs(corrected_dist - dist)="<<fabs(corrected_dist - dist)<<"  cdchit->dist="<<cdchit->dist<<"  dist="<<dist<<endl;
			if(!match)continue;
			
			// Check if this is already in the map
			map<const DCDCTrackHit*, const DMCTrackHit*>::iterator iter = cdclink.find(cdchit);
			if(iter!=cdclink.end()){
				// If more than one hit occurs in this straw then keep the
				// one closest to corrected_dist
				const DMCTrackHit *mchit_old = cdclink[cdchit];
				DVector3 pos_truth(mchit_old->r*cos(mchit_old->phi), mchit_old->r*sin(mchit_old->phi), mchit_old->z);
				DVector3 A = wire->udir.Cross(pos_truth - wire->origin);
				double dist_old = A.Mag();
				if(fabs(corrected_dist - dist)<fabs(corrected_dist - dist_old)){
					cdclink[cdchit] = mchit;
				}
			}else{
				cdclink[cdchit] = mchit;
			}
		}
	}
}



