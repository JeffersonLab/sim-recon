// $Id: DEventProcessor_mcthrown_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_mcthrown_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_mcthrown_hists.h"

#include "DApplication.h"
#include "DEventLoop.h"
#include "DMCThrown.h"

static TFile **tfilePtr = NULL;

// Routine used to create our DEventProcessor
extern "C"{
void InitProcessors(DApplication *app){
	app->AddProcessor(new DEventProcessor_mcthrown_hists());
}

void SetTFilePtrAddress(TFile **h){
	tfilePtr = h;
}
} // "C"


//------------------
// DEventProcessor_mcthrown_hists
//------------------
DEventProcessor_mcthrown_hists::DEventProcessor_mcthrown_hists()
{	
}

//------------------
// ~DEventProcessor_mcthrown_hists
//------------------
DEventProcessor_mcthrown_hists::~DEventProcessor_mcthrown_hists()
{
}

//------------------
// init
//------------------
derror_t DEventProcessor_mcthrown_hists::init(void)
{
	// open ROOT file (if needed)
	ROOTfile = NULL;
	if(tfilePtr == NULL)tfilePtr = &ROOTfile;
	if(*tfilePtr == NULL){
		*tfilePtr = ROOTfile = new TFile("mcthrown_hists.root","RECREATE","Produced by hd_ana");
		cout<<"Opened ROOT file \"mcthrown_hists.root\""<<endl;
	}else{
		(*tfilePtr)->cd();
	}

	// Create THROWN directory
	TDirectory *dir = new TDirectory("THROWN","THROWN");
	dir->cd();

	// Create histograms
	pmom	= new TH1F("pmom","Thrown momentum in GeV/c",1200, 0.0, 12.0);
	theta	= new TH1F("theta","Thrown theta in radians",200, 0.0, M_PI);
	phi	= new TH1F("phi","Thrown phi in radians",200, 0.0, 2.0*M_PI);
	energy	= new TH1F("energy","Thrown energy in GeV",1200, 0.0, 12.0);

	Nparticles_per_event	= new TH1F("Nparticles_per_event","Number of thrown particles per event",21, -0.5, 20.5);
	particle_type	= new TH1F("particle_type","GEANT3 particle type of thrown particles",101, -0.5, 100.5);
	
	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
derror_t DEventProcessor_mcthrown_hists::evnt(DEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	
	// Loop over thrown tracks
	Nparticles_per_event->Fill(mcthrowns.size());
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];

		pmom->Fill(mcthrown->p);
		theta->Fill(mcthrown->theta);
		phi->Fill(mcthrown->phi);
		energy->Fill(mcthrown->E);
		vertex->Fill(mcthrown->x, mcthrown->y, mcthrown->z);
		
		particle_type->Fill(mcthrown->type);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
derror_t DEventProcessor_mcthrown_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
derror_t DEventProcessor_mcthrown_hists::fini(void)
{

	if(ROOTfile){
		ROOTfile->Write();
		delete ROOTfile;
		cout<<endl<<"Closed ROOT file"<<endl;
	}

	return NOERROR;
}
