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

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
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
	if(ROOTfile != NULL) ROOTfile->cd();
	
	// Create THROWN directory
	TDirectory *dir = new TDirectory("CDC","CDC");
	dir->cd();

	// Create histograms
	dE	= new TH1F("dE","Energy loss in staw tube in keV",1000, 0.0, 100.0);
	
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
	loop->Get(cdchits);
	
	// Loop over thrown tracks
	for(unsigned int i=0;i<cdchits.size();i++){
		const DCDCHit *cdchit = cdchits[i];

		dE->Fill(cdchit->dE*1.0E6);
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
