// $Id: DEventProcessor_mcthrown_hists.cc 1816 2006-06-06 14:38:18Z davidl $
//
//    File: DEventProcessor_mcthrown_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include "DEventProcessor_mcthrown_hists.h"
#include "TRACKING/DMCThrown.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our JEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_mcthrown_hists());
}
}

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
jerror_t DEventProcessor_mcthrown_hists::init(void)
{
	// open ROOT file (if needed)
	//if(ROOTfile != NULL) ROOTfile->cd();

	// Create THROWN directory
	TDirectory *dir = new TDirectory("THROWN","THROWN");
	dir->cd();

	// Create histograms
	pmom	= new TH1F("pmom","Thrown momentum in GeV/c",1200, 0.0, 12.0);
	theta	= new TH1F("theta","Thrown theta in degrees",1000, 0.0, 180.0);
	phi	= new TH1F("phi","Thrown phi in radians",200, 0.0, 2.0*M_PI);
	energy	= new TH1F("energy","Thrown energy in GeV",1200, 0.0, 12.0);
	pmom_vs_theta = new TH2F("pmom_vs_theta","Thrown momentum vs. #theta",1000, 0.0, 180.0, 300, 0.0, 10.0);
	pmom_vs_theta_pip = new TH2F("pmom_vs_theta_pip","Thrown momentum vs. #theta for #pi^{+}",1000, 0.0, 180.0, 300, 0.0, 10.0);
	pmom_vs_theta_pim = new TH2F("pmom_vs_theta_pim","Thrown momentum vs. #theta for #pi^{-}",1000, 0.0, 180.0, 300, 0.0, 10.0);
	pmom_vs_theta_proton = new TH2F("pmom_vs_theta_proton","Thrown momentum vs. #theta for protons",1000, 0.0, 180.0, 300, 0.0, 10.0);
	pmom_vs_theta_gamma = new TH2F("pmom_vs_theta_gamma","Thrown momentum vs. theta for gammas",1000, 0.0, 180.0, 300, 0.0, 10.0);

	vertex = new TH3F("vertex", "Position of vertex from which thrown particles were thrown", 50, -5.0, 5.0, 50, -5.0, 5.0, 150, 0.0, 150.0);

	Nparticles_per_event	= new TH1F("Nparticles_per_event","Number of thrown particles per event",21, -0.5, 20.5);
	particle_type	= new TH1F("particle_type","GEANT3 particle type of thrown particles",101, -0.5, 100.5);
	
	// Go back up to the parent directory
	dir->cd("../");
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_mcthrown_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCThrown*> mcthrowns;
	loop->Get(mcthrowns);
	
	// Loop over thrown tracks
	Nparticles_per_event->Fill(mcthrowns.size());
	for(unsigned int i=0;i<mcthrowns.size();i++){
		const DMCThrown *mcthrown = mcthrowns[i];

		pmom->Fill(mcthrown->p);
		theta->Fill(mcthrown->theta*57.3);
		phi->Fill(mcthrown->phi);
		energy->Fill(mcthrown->E);
		pmom_vs_theta->Fill(mcthrown->theta*57.3, mcthrown->p);
		vertex->Fill(mcthrown->x, mcthrown->y, mcthrown->z);
		
		particle_type->Fill(mcthrown->type);
		
		switch(mcthrown->type){
			case 1:
				pmom_vs_theta_gamma->Fill(mcthrown->theta*57.3, mcthrown->p);
				break;
			case 8:
				pmom_vs_theta_pip->Fill(mcthrown->theta*57.3, mcthrown->p);
				break;
			case 9:
				pmom_vs_theta_pim->Fill(mcthrown->theta*57.3, mcthrown->p);
				break;
			case 14:
				pmom_vs_theta_proton->Fill(mcthrown->theta*57.3, mcthrown->p);
				break;
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_mcthrown_hists::erun(void)
{

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_mcthrown_hists::fini(void)
{
	return NOERROR;
}
