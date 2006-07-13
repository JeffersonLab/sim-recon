// $Id: DEventProcessor_acceptance_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_acceptance_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <map>
using namespace std;

#include <TThread.h>

#include <JANA/JEventLoop.h>

#include "DEventProcessor_acceptance_hists.h"

#include "DANA/DApplication.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackHit.h"
#include "TRACKING/DTrackHit_factory_MC.h"

#define MIN_CDC_HITS 8
#define MIN_FDC_HITS 8
#define MIN_CDC_FDC_HITS 8
#define MIN_BCAL_HITS 4
#define MIN_FCAL_HITS 4
#define MIN_TOF_HITS 1
#define MIN_UPV_HITS 4

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
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
	if(ROOTfile != NULL) ROOTfile->cd();

	// Create ACCEPTANCE directory
	TDirectory *dir = new TDirectory("ACCEPTANCE","ACCEPTANCE");
	dir->cd();

	// Create histograms
	int N_p_bins = 100;
	float p_lo = 0.0;
	float p_hi = 12.0;
	int N_theta_bins = 100;
	float theta_lo = 0.0;
	float theta_hi = M_PI;
	CDC	= new TH2F("CDC","CDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	FDC	= new TH2F("FDC","FDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	CDC_FDC	= new TH2F("CDC_FDC","CDC_FDC acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	BCAL	= new TH2F("BCAL","BCAL acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	FCAL	= new TH2F("FCAL","FCAL acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	TOF	= new TH2F("TOF","TOF acceptance",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	thrown_charged	= new TH2F("thrown_charged","thrown charged particles",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	thrown_photon	= new TH2F("thrown_photon","thrown photons",N_p_bins, p_lo, p_hi, N_theta_bins, theta_lo, theta_hi);
	
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
	vector<const DTrackHit*> trackhits;
	loop->Get(mcthrowns);
	JFactory_base *fac = loop->Get(trackhits, "MC");
	DTrackHit_factory_MC *factory = dynamic_cast<DTrackHit_factory_MC*>(fac);
	
	// Loop over thrown tracks
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];
		
		if(mcthrown->q != 0.0){
			thrown_charged->Fill(mcthrown->p, mcthrown->theta);
		}else if(mcthrown->type == 1){
			thrown_photon->Fill(mcthrown->p, mcthrown->theta);
		}
	}
	
	// Make sure we got a valid pointer to the trackhits MC factory.
	if(!factory){
		static int warnings=0;
		if(warnings<10){
			cout<<__FILE__<<":"<<__LINE__<<" Unable to get point to DFactory_DTrackHit_MC object!"<<endl;	
			warnings++;
		}
		if(warnings == 10)cout<<__FILE__<<":"<<__LINE__<<" LAST WARNING!!"<<endl;
		return NOERROR;
	}

	// Loop over track hits
	map<int,int> cdchits;
	map<int,int> fdchits;
	map<int,int> bcalhits;
	const vector<int>& tracknumber = factory->GetTrackNumbers();
	const vector<bool>& primaryflag = factory->GetPrimaryFlags();
	for(unsigned int i=0;i<trackhits.size();i++){
		const DTrackHit *trackhit = trackhits[i];

		switch(trackhit->system){
			case SYS_CDC:
				if(primaryflag[i])cdchits[tracknumber[i]]++;
				break;
			case SYS_FDC:
				if(primaryflag[i])fdchits[tracknumber[i]]++;
				break;
			case SYS_BCAL:
				bcalhits[tracknumber[i]]++;
				break;
			default:
				break;
		}
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
		
		if(iter->second >= MIN_CDC_HITS)CDC->Fill(mcthrown->p, mcthrown->theta);
	}

	// Loop over tracks in the FDC
	for(iter=fdchits.begin(); iter!=fdchits.end(); iter++){

		// Find thrown parameters for this track (if any)
		if(iter->first<=0 || iter->first>(int)mcthrowns.size())continue;
		const DMCThrown *mcthrown = mcthrowns[iter->first-1];
		
		if(iter->second >= MIN_FDC_HITS)FDC->Fill(mcthrown->p, mcthrown->theta);
		
		// Fill CDC + FDC histo
		if(cdchits.find(iter->first) != cdchits.end()){
			int cdc_fdc_hits = cdchits.find(iter->first)->second + iter->second;
			if(cdc_fdc_hits >= MIN_CDC_FDC_HITS)
				CDC_FDC->Fill(mcthrown->p, mcthrown->theta);
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
