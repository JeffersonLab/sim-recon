// $Id: $
//
//    File: DEventProcessor_cdc_hists.cc
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include <JANA/JEventLoop.h>

#include "DEventProcessor_cdc_hists.h"

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DTrack.h>
#include <TRACKING/DMCTrackHit.h>
#include <CDC/DCDCHit.h>
#include <DVector2.h>


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
	cdc_ptr = &cdc;
	cdchit_ptr = &cdchit;
	
	pthread_mutex_init(&mutex, NULL);
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

	// Create Tree
	cdctree = new TTree("cdc","CDC Truth points");
	cdchittree = new TTree("cdchit","CDC Hits");
	cdcbranch = cdctree->Branch("T","CDC_branch",&cdc_ptr);
	cdchitbranch = cdchittree->Branch("H","CDChit_branch",&cdchit_ptr);

	// Go back up to the parent directory
	dir->cd("../");
	
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

//------------------
// evnt
//------------------
jerror_t DEventProcessor_cdc_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DCDCHit*> cdchits;
	loop->Get(mctrackhits);
	loop->Get(cdchits);
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Loop over all truth hits, ignoring all but FDC hits
	for(unsigned int i=0; i<mctrackhits.size(); i++){
		const DMCTrackHit *mctrackhit = mctrackhits[i];
		if(mctrackhit->system != SYS_CDC)continue;
		
		double r = mctrackhit->r;
		double phi = mctrackhit->phi;
		double x = r*cos(phi);
		double y = r*sin(phi);
		cdc.pos_truth.SetXYZ(x, y, mctrackhit->z);
		
		cdctree->Fill();
	}

	// Loop over all real hits
	for(unsigned int i=0; i<cdchits.size(); i++){
		const DCDCHit *hit = cdchits[i];
		
		cdchit.ring		= hit->ring;
		cdchit.straw	= hit->straw;
		cdchit.dE		= hit->dE;
		cdchit.t			= hit->t;
		
		cdchittree->Fill();
	}

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}
