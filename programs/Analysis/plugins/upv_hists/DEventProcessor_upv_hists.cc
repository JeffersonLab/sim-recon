// $Id$
//
//    File: DEventProcessor_upv_hists.cc
// Created: Mon Apr  3 11:38:03 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_upv_hists.h"

#include <TLorentzVector.h>

#include "DANA/DApplication.h"
#include "UPV/DUPVTruthHit_factory.h"
#include "TRACKING/DMCThrown.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_upv_hists());
}
} // "C"


#define FCAL_Z_OFFSET 640.0-65.0 // I don't know what this value is ???
//#define FCAL_Z_OFFSET 170.0 // I don't know what this value is ???
#define PI_ZERO_MASS 0.13497

//------------------
// init
//------------------
jerror_t DEventProcessor_upv_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();
	
	xy_shower = new TH2F("upv_xy_shower","xy_shower",100, -100.0, 100., 100 , -100.0, 100.0);
	z_shower = new TH1F("upv_z_shower","z_shower",1000, -120.0, 20.);
	E_shower = new TH1F("upv_E_shower","E_shower", 4000, 0.0, 0.010);
	Etotal = new TH1F("upv_Etotal","Etotal", 4000, 0.0, 1.000);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_upv_hists::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_upv_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DUPVTruthHit*> showers;
	DUPVTruthHit_factory *fac = dynamic_cast<DUPVTruthHit_factory*>(loop->Get(showers));
	
	LockState();
	
	// Single shower params
	for(unsigned int i=0; i<showers.size(); i++){
		const DUPVTruthHit *s = showers[i];
		xy_shower->Fill(s->x, s->y);
		z_shower->Fill(s->z);
		E_shower->Fill(s->E);
	}
	
	if(fac)Etotal->Fill(fac->GetETotal());

	UnlockState();	

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_upv_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_upv_hists::fini(void)
{
	return NOERROR;
}

