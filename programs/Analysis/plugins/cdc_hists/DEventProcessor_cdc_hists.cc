// $Id: $
//
//    File: DEventProcessor_cdc_hists.cc
//

#include <iostream>
using namespace std;

#include <TThread.h>

#include <JANA/JEventLoop.h>

#include "DEventProcessor_cdc_hists.h"

#include "DANA/DApplication.h"
#include "CDC/DCDCHit.h"
#include "TRACKING/DMCThrown.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_cdc_hists());
}
} // "C"


//------------------
// DEventProcessor_cdc_hists
//------------------
DEventProcessor_cdc_hists::DEventProcessor_cdc_hists()
{	
}

//------------------
// ~DEventProcessor_cdc_hists
//------------------
DEventProcessor_cdc_hists::~DEventProcessor_cdc_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_cdc_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();
	
	// Create THROWN directory
	TDirectory *dir = new TDirectory("CDC","CDC");
	dir->cd();

	// Create histograms
	dE	= new TH1F("dE","Energy loss in staw tube in keV",1000, 0.0, 100.0);
	cdc_layer1_theta_vs_p = new TH2F("cdc_layer1_theta_vs_p","#theta vs. momentum for tracks hitting CDC layer 1",25, 0.1, 5.1, 200, 0.0, 40.0);
	cdc_layer23_theta_vs_p = new TH2F("cdc_layer23_theta_vs_p","#theta vs. momentum for tracks hitting CDC layer 23",25, 0.1, 5.1, 200, 0.0, 40.0);
	
	cdc_nhits_per_event = new TH1F("cdc_nhits_per_event", "CDC hits per event", 501, 0.0, 500.0);

	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_cdc_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DCDCHit*> cdchits;
	vector<const DMCThrown*> mcthrowns;
	loop->Get(cdchits);
	loop->Get(mcthrowns);
	
	cdc_nhits_per_event->Fill((double)cdchits.size());
	
	// Loop over CDC hits
	bool ishit1 = false;
	bool ishit23 = false;
	for(unsigned int i=0;i<cdchits.size();i++){
		const DCDCHit *cdchit = cdchits[i];

		dE->Fill(cdchit->dE*1.0E6);
		
		if(cdchit->ring == 1)ishit1 = true;
		if(cdchit->ring == 23)ishit23 = true;
	}
	
	if(mcthrowns.size()==1){
		double p = mcthrowns[0]->p;
		double theta = mcthrowns[0]->theta*57.3;
		if(ishit1)cdc_layer1_theta_vs_p->Fill(p,theta);
		if(ishit23)cdc_layer23_theta_vs_p->Fill(p,theta);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_cdc_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_cdc_hists::fini(void)
{
	return NOERROR;
}
