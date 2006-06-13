// $Id: $
//
//    File: DEventProcessor_cdc_hists.cc
//

#include <iostream>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_cdc_hists.h"

#include "DApplication.h"
#include "DEventLoop.h"
#include "DCDCHit.h"

static TFile **tfilePtr = NULL;

// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
	app->AddProcessor(new DEventProcessor_cdc_hists());
}

void SetTFilePtrAddress(TFile **h){
	tfilePtr = h;
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
derror_t DEventProcessor_cdc_hists::init(void)
{
	// open ROOT file (if needed)
	ROOTfile = NULL;
	if(tfilePtr == NULL)tfilePtr = &ROOTfile;
	if(*tfilePtr == NULL){
		*tfilePtr = ROOTfile = new TFile("cdc_hists.root","RECREATE","Produced by hd_ana");
		cout<<"Opened ROOT file \"cdc_hists.root\""<<endl;
	}else{
		(*tfilePtr)->cd();
	}

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
derror_t DEventProcessor_cdc_hists::evnt(DEventLoop *loop, int eventnumber)
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
derror_t DEventProcessor_cdc_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
derror_t DEventProcessor_cdc_hists::fini(void)
{

	if(ROOTfile){
		ROOTfile->Write();
		delete ROOTfile;
		cout<<endl<<"Closed ROOT file"<<endl;
	}

	return NOERROR;
}
