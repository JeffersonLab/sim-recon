// $Id: DEventProcessor_fdc_hists.cc 2774 2007-07-19 15:59:02Z davidl $
//
//    File: DEventProcessor_fdc_hists.cc
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <iostream>
#include <cmath>
using namespace std;

#include <TThread.h>

#include "DEventProcessor_fdc_hists.h"

#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include <DANA/DApplication.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrackHit.h>
#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCHit.h>
#include <DVector2.h>

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
//extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new DEventProcessor_fdc_hists());
}
} // "C"


//------------------
// DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::DEventProcessor_fdc_hists()
{
	fdc_ptr = &fdc;
	fdchit_ptr = &fdchit;
	
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~DEventProcessor_fdc_hists
//------------------
DEventProcessor_fdc_hists::~DEventProcessor_fdc_hists()
{
}

//------------------
// init
//------------------
jerror_t DEventProcessor_fdc_hists::init(void)
{
	// Create TRACKING directory
	TDirectory *dir = new TDirectoryFile("FDC","FDC");
	dir->cd();

	// Create Tree
	fdctree = new TTree("fdc","FDC Truth points");
	fdchittree = new TTree("fdchit","FDC Hits");
	fdcbranch = fdctree->Branch("T","FDC_branch",&fdc_ptr);
	fdchitbranch = fdchittree->Branch("H","FDChit_branch",&fdchit_ptr);

	dir->cd("../");
	

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_fdc_hists::brun(JEventLoop *loop, int runnumber)
{	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_fdc_hists::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_fdc_hists::fini(void)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_fdc_hists::evnt(JEventLoop *loop, int eventnumber)
{
	vector<const DMCTrackHit*> mctrackhits;
	vector<const DFDCHit*> fdchits;
	loop->Get(mctrackhits);
	loop->Get(fdchits);
	
	// Lock mutex
	pthread_mutex_lock(&mutex);
	
	// Loop over all truth hits, ignoring all but FDC hits
	for(unsigned int i=0; i<mctrackhits.size(); i++){
		const DMCTrackHit *mctrackhit = mctrackhits[i];
		if(mctrackhit->system != SYS_FDC)continue;
		
		double r = mctrackhit->r;
		double phi = mctrackhit->phi;
		double x = r*cos(phi);
		double y = r*sin(phi);
		fdc.pos_truth.SetXYZ(x, y, mctrackhit->z);
		
		fdctree->Fill();
	}
#if 1
	// Loop over all real hits
	for(unsigned int i=0; i<fdchits.size(); i++){
		const DFDCHit *hit = fdchits[i];
		
		fdchit.layer	= hit->layer;
		fdchit.module	= hit->module;
		fdchit.element	= hit->element;
		fdchit.plane	= hit->plane;
		fdchit.gPlane	= hit->gPlane;
		fdchit.gLayer	= hit->gLayer;
		fdchit.q			= hit->q;
		fdchit.t			= hit->t;
		fdchit.r			= hit->r;
		fdchit.type		= hit->type;
		
		fdchittree->Fill();
	}
#endif

	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

