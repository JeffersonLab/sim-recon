// $Id: DEventProcessor_acceptance_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_acceptance_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <map>
#include <cmath>
using namespace std;

#include <TThread.h>

#include <JANA/JEventLoop.h>

#include "DEventProcessor_acceptance_hists.h"

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCGeometry.h>
#include <CDC/DCDCTrackHit.h>

#define MIN_CDC_HITS 8
#define MIN_FDC_HITS 8
#define MIN_CDC_FDC_HITS 8
#define MIN_BCAL_HITS 4
#define MIN_FCAL_HITS 4
#define MIN_TOF_HITS 1
#define MIN_UPV_HITS 4

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_acceptance_hists());
}
} // "C"


//------------------
// DEventProcessor_acceptance_hists
//------------------
DEventProcessor_acceptance_hists::DEventProcessor_acceptance_hists()
{	
}

//------------------
// ~DEventProcessor_acceptance_hists
//------------------
DEventProcessor_acceptance_hists::~DEventProcessor_acceptance_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_acceptance_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();

	// Create ACCEPTANCE directory
	TDirectory *dir = new TDirectoryFile("ACCEPTANCE","ACCEPTANCE");
	dir->cd();

	// Create histograms
	int N_p_bins = 100;
	float p_lo = 0.0;
	float p_hi = 12.0;
	int N_theta_bins = 400;
	float theta_lo = 0.0;
	float theta_hi = M_PI*57.3;
	CDC	= new TH2F("CDC","CDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	FDC	= new TH2F("FDC","FDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	CDC_FDC	= new TH2F("CDC_FDC","CDC_FDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	BCAL	= new TH2F("BCAL","BCAL acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	FCAL	= new TH2F("FCAL","FCAL acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	TOF	= new TH2F("TOF","TOF acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	thrown_charged	= new TH2F("thrown_charged","thrown charged particles",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	thrown_photon	= new TH2F("thrown_photon","thrown photons",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	
	FDC_anode_hits_per_event = new TH1D("FDC_anode_hits_per_event","FDC anode hits/event", 201, -0.5, 200.5);
	FDC_anode_hits_per_layer = new TH1D("FDC_anode_hits_per_layer","FDC anode hits/layer", 24, 0.5, 24.5);
	FDC_anode_hits_per_wire = new TH1D("FDC_anode_hits_per_wire","FDC anode hits/wire", 96, 0.5, 96.5);

	CDC_nhits_vs_pthrown = new TH1D("CDC_nhits_vs_pthrown","Number of CDC hits per event vs. thrown momentum", 40, 0.0, 9.0);
	FDC_nhits_vs_pthrown = new TH1D("FDC_nhits_vs_pthrown","Number of FDC anode hits per event vs. thrown momentum", 40, 0.0, 9.0);
	pthrown = new TH1D("pthrown","thrown momentum", 40, 0.0, 9.0);
	
	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_acceptance_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	vector<const DMCTrackHit*> mctrackhits;
	loop->Get(mcthrowns);
	loop->Get(mctrackhits);	
	
	// Count FDC anode hits
	double Nfdc_anode = 0.0;
	vector<const DFDCHit*> fdchits;
	loop->Get(fdchits);
	for(unsigned int i=0; i<fdchits.size(); i++){
		if(fdchits[i]->type==0){
			Nfdc_anode+=1.0;
			
			FDC_anode_hits_per_layer->Fill(DFDCGeometry::gLayer(fdchits[i]));
			FDC_anode_hits_per_wire->Fill(fdchits[i]->element);
		}
	}
	FDC_anode_hits_per_event->Fill(Nfdc_anode);
	
	// Count CDC hits
	double Ncdc_anode = 0.0;
	vector<const DCDCTrackHit*> cdctrackhits;
	loop->Get(cdctrackhits);
	Ncdc_anode = (double)cdctrackhits.size();
	
	
	// Loop over thrown tracks
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		if(mcthrown->charge() != 0.0){
			thrown_charged->Fill(mcthrown->momentum().Mag(), mcthrown->momentum().Theta()*57.3);
		}else if(mcthrown->type == 1){
			thrown_photon->Fill(mcthrown->momentum().Mag(), mcthrown->momentum().Theta()*57.3);
		}
	}
	
	// Loop over track hits
	map<int,int> cdchits;
	map<int,int> fdchitmap;
	map<int,int> bcalhits;
	for(unsigned int i=0;i<mctrackhits.size();i++){
		const DMCTrackHit *mctrackhit = mctrackhits[i];

		switch(mctrackhit->system){
			case SYS_CDC:
				if(mctrackhit->primary)cdchits[mctrackhit->track]++;
				break;
			case SYS_FDC:
				if(mctrackhit->primary)fdchitmap[mctrackhit->track]++;
				break;
			case SYS_BCAL:
				bcalhits[mctrackhit->track]++;
				break;
			default:
				break;
		}
	}
	
	// Simple 1-D histos for CDC and FDC as a function of thrown momentum
	if(mcthrowns.size()==1){
		double p = mcthrowns[0]->momentum().Mag();
		CDC_nhits_vs_pthrown->Fill(p, Ncdc_anode);
		FDC_nhits_vs_pthrown->Fill(p, Nfdc_anode);
		pthrown->Fill(p);
	}

	
	// NOTE: In the sections that follow we have to assume that
	// the track number of the hit corresponds to the position of
	// the DMCThrown object in the list of DMCThrown objects.
	// This should be a good assumption, but I don't know that
	// it is (or always will be) guaranteed.
	
	// Loop over tracks in the CDC
	map<int,int>::iterator iter;
	for(iter=cdchits.begin(); iter!=cdchits.end(); iter++){

		// Find thrown parameters for this track (if any)
		if(iter->first<=0 || iter->first>(int)mcthrowns.size())continue;
		const DMCThrown *mcthrown = mcthrowns[iter->first-1];
		
		if(iter->second >= MIN_CDC_HITS)CDC->Fill(mcthrown->momentum().Mag(), mcthrown->momentum().Theta()*57.3);
	}

	// Loop over tracks in the FDC
	for(iter=fdchitmap.begin(); iter!=fdchitmap.end(); iter++){

		// Find thrown parameters for this track (if any)
		if(iter->first<=0 || iter->first>(int)mcthrowns.size())continue;
		const DMCThrown *mcthrown = mcthrowns[iter->first-1];
		
		if(iter->second >= MIN_FDC_HITS)FDC->Fill(mcthrown->momentum().Mag(), mcthrown->momentum().Theta()*57.3);
		
		// Fill CDC + FDC histo
		if(cdchits.find(iter->first) != cdchits.end()){
			int cdc_fdc_hits = cdchits.find(iter->first)->second + iter->second;
			if(cdc_fdc_hits >= MIN_CDC_FDC_HITS)
				CDC_FDC->Fill(mcthrown->momentum().Mag(), mcthrown->momentum().Theta()*57.3);
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_acceptance_hists::erun(void)
{
	CDC->Divide(thrown_charged);
	FDC->Divide(thrown_charged);
	CDC_FDC->Divide(thrown_charged);

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_acceptance_hists::fini(void)
{

	return NOERROR;
}
