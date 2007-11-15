// $Id: DEventProcessor_trackeff_hists.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_trackeff_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_trackeff_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <FDC/DFDCGeometry.h>

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_trackeff_hists());
}
} // "C"


//------------------
// DEventProcessor_trackeff_hists
//------------------
DEventProcessor_trackeff_hists::DEventProcessor_trackeff_hists()
{
	leaf_ptr = &leaf;
	
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&rt_mutex, NULL);
}

//------------------
// ~DEventProcessor_trackeff_hists
//------------------
DEventProcessor_trackeff_hists::~DEventProcessor_trackeff_hists()
{
	delete ref;
}

//------------------
// init
//------------------
jerror_t DEventProcessor_trackeff_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = new TDirectory("TRACKING","TRACKING");
	dir->cd();

	// Create Tree
	trkeff = new TTree("trkeff","Tracking Efficiency");
	trkeff->Branch("F","TrkEff_Leaf",&leaf_ptr);

	dir->cd("../");
	
	MAX_HIT_DIST_CDC = 1.0; // cm
	MAX_HIT_DIST_FDC = 5.0; // cm

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_trackeff_hists::brun(JEventLoop *loop, int runnumber)
{
	DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
	bfield = dapp->GetBfield(); // temporary until new geometry scheme is worked out
	ref = new DReferenceTrajectory(bfield);
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_trackeff_hists::evnt(JEventLoop *loop, int eventnumber)
{
	CDChitv cdctrackhits;
	FDChitv fdctrackhits;
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DMCThrown*> mcthrowns;
	
	loop->Get(cdctrackhits);
	loop->Get(fdctrackhits);
	loop->Get(trackcandidates);
	loop->Get(mcthrowns);
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Get hit list for all candidates
	vector<CDChitv> cdc_candidate_hits;
	vector<FDChitv> fdc_candidate_hits;
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		CDChitv cdc_outhits;
		GetCDCHits(trackcandidates[i], cdctrackhits, cdc_outhits);
		cdc_candidate_hits.push_back(cdc_outhits);

		FDChitv fdc_outhits;
		GetFDCHits(trackcandidates[i], fdctrackhits, fdc_outhits);
		fdc_candidate_hits.push_back(fdc_outhits);
	}
	
	int Nnon_noise_cdcs = (int)cdctrackhits.size() - 100; // hardwire approximate value for now
	if(Nnon_noise_cdcs<0)Nnon_noise_cdcs=0;
	leaf.status = ((double)Nnon_noise_cdcs/25.0 - mcthrowns.size())<1.5 ? 0:-1;

	// Get hit list for all throwns
	for(unsigned int i=0; i<mcthrowns.size(); i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		// if this isn't a charged track, then skip it
		if(fabs(mcthrowns[i]->charge())==0.0)continue;

		CDChitv cdc_thrownhits;
		FDChitv fdc_thrownhits;
		GetCDCHits(mcthrowns[i], cdctrackhits, cdc_thrownhits);
		GetFDCHits(mcthrowns[i], fdctrackhits, fdc_thrownhits);
		
		if(cdc_thrownhits.size()<5 && fdc_thrownhits.size()<5) leaf.status--;
		
		leaf.pthrown = mcthrown->momentum();
		leaf.z_thrown = mcthrown->position().Z();
		leaf.ncdc_hits_thrown = cdc_thrownhits.size();
		leaf.nfdc_hits_thrown = fdc_thrownhits.size();
		leaf.ncdc_hits = cdctrackhits.size();
		leaf.nfdc_hits = GetNFDCWireHits(fdctrackhits);
		
		// Look for a candidate track that matches this thrown one
		CDChitv cdc_matched_hits;
		FDChitv fdc_matched_hits;
		unsigned int icdc_can = FindMatch(cdc_thrownhits, cdc_candidate_hits, cdc_matched_hits);
		unsigned int ifdc_can = FindMatch(fdc_thrownhits, fdc_candidate_hits, fdc_matched_hits);
		
		// Initialize values for when no matching candidate is found
		leaf.pcan.SetXYZ(0.0, 0.0, 0.0);
		leaf.z_can = -1000.0;
		leaf.ncdc_hits_can = 0;
		leaf.nfdc_hits_can = 0;
		leaf.ncdc_hits_thrown_and_can = 0;
		leaf.nfdc_hits_thrown_and_can = 0;
		leaf.cdc_chisq_can = 1.0E6;
		leaf.fdc_chisq_can = 1.0E6;
		
		// CDC
		const DTrackCandidate *cdc_can = NULL;
		if(icdc_can>=0 && icdc_can<trackcandidates.size()){
			cdc_can = trackcandidates[icdc_can];
			
			leaf.ncdc_hits_can = cdc_candidate_hits[icdc_can].size();
			leaf.ncdc_hits_thrown_and_can = cdc_matched_hits.size();
			leaf.cdc_chisq_can = 0.0;
		}
		
		// FDC
		const DTrackCandidate *fdc_can = NULL;
		if(ifdc_can>=0 && ifdc_can<trackcandidates.size()){
			fdc_can = trackcandidates[ifdc_can];
			
			leaf.nfdc_hits_can = fdc_candidate_hits[ifdc_can].size();
			leaf.nfdc_hits_thrown_and_can = fdc_matched_hits.size();
			leaf.fdc_chisq_can = 0.0;
		}
		
		// Figure out which candidate (if any) I should match this with
		const DTrackCandidate* can = NULL;
		if(cdc_can!=NULL && fdc_can!=NULL){
			can = cdc_matched_hits.size()>fdc_matched_hits.size() ? cdc_can:fdc_can;
		}else{
			can = cdc_can!=NULL ? cdc_can:fdc_can;
		}
		
		if(can!=NULL){
			leaf.pcan = can->momentum();
			leaf.z_can = can->position().Z();
		}
		
		trkeff->Fill();
	}
	
	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// GetCDCHits
//------------------
void DEventProcessor_trackeff_hists::GetCDCHits(const DKinematicData *p, CDChitv &inhits, CDChitv &outhits)
{
	// In case we run with multiple threads
	pthread_mutex_lock(&rt_mutex);
	
	// Re-swim the reference trajectory using the parameters in p
	ref->Swim(p->position(), p->momentum(), p->charge());
	
	// Loop over hits in the "inhits" vector and find the DOCA of the R.T.
	// for each.
	outhits.clear();
	for(unsigned int i=0; i<inhits.size(); i++){
		double doca = ref->DistToRT(inhits[i]->wire);
		if(doca < MAX_HIT_DIST_CDC)outhits.push_back(inhits[i]);
	}
	
	pthread_mutex_unlock(&rt_mutex);
}

//------------------
// GetFDCHits
//------------------
void DEventProcessor_trackeff_hists::GetFDCHits(const DKinematicData *p, FDChitv &inhits, FDChitv &outhits)
{
	// In case we run with multiple threads
	pthread_mutex_lock(&rt_mutex);
	
	// Re-swim the reference trajectory using the parameters in p
	ref->Swim(p->position(), p->momentum(), p->charge());
	
	// Loop over hits in the "inhits" vector and find the DOCA of the R.T.
	// for each.
	outhits.clear();
	for(unsigned int i=0; i<inhits.size(); i++){
		const DFDCHit* hit = inhits[i];
		if(hit->type != 0)continue; // filter out cathode hits
		const DFDCWire *wire = DFDCGeometry::GetDFDCWire(hit->gLayer, hit->element);
		if(!wire)continue;
		double doca = ref->DistToRT(wire);
		if(doca < MAX_HIT_DIST_FDC)outhits.push_back(hit);
	}
	
	pthread_mutex_unlock(&rt_mutex);
}

//------------------
// GetNFDCWireHits
//------------------
unsigned int DEventProcessor_trackeff_hists::GetNFDCWireHits(FDChitv &inhits)
{
	unsigned int N=0;
	for(unsigned int i=0; i<inhits.size(); i++){
		if(inhits[i]->type==0)N++;
	}
	return N;
}

//------------------
// FindMatch
//------------------
unsigned int DEventProcessor_trackeff_hists::FindMatch(CDChitv &thrownhits, vector<CDChitv> &candidate_hits, CDChitv &matched_hits)
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
// FindMatch
//------------------
unsigned int DEventProcessor_trackeff_hists::FindMatch(FDChitv &thrownhits, vector<FDChitv> &candidate_hits, FDChitv &matched_hits)
{
	// Loop through all of the candidate hits and look for the one that best
	// matches this thrown's hits. The algorithm simply looks for the candidate
	// that has the most hits in common with the thrown.
	//
	// If the best match shares less than half its hits with the thrown,
	// then it is considered not to match and -1 is returned.

	unsigned int ibest=(unsigned int)-1;
	FDChitv my_matched;
	matched_hits.clear();
	for(unsigned int i=0; i<candidate_hits.size(); i++){
		FDChitv &canhits = candidate_hits[i];
		
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
//_DBG_<<"matched_hits.size()="<<matched_hits.size()<<"  thrownhits.size()="<<thrownhits.size()<<endl;
	if(matched_hits.size() >= (thrownhits.size()/2))return ibest;
	
	return (unsigned int)-1;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_trackeff_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_trackeff_hists::fini(void)
{
	return NOERROR;
}
