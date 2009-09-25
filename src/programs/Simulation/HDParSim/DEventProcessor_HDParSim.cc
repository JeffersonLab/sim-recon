// $Id$
//
//    File: DEventProcessor_HDParSim.cc
// Created: Tue Feb  3 08:39:28 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include "DEventProcessor_HDParSim.h"


//------------------
// DEventProcessor_HDParSim  (constructor)
//------------------
DEventProcessor_HDParSim::DEventProcessor_HDParSim(const char *fname)
{
	pthread_mutex_init(&root_mutex, NULL);
	rootFile = NULL;
	rootFileName = fname;

	// Create a DTrackingResolution object. Specifically,
	// we create a DTrackingResolutionGEANT object so as
	// to use the resolutions derived from the GEANT-based
	// simulations.
	res = new DTrackingResolutionGEANT();
}

//------------------
// init
//------------------
jerror_t DEventProcessor_HDParSim::init(void)
{
	pthread_mutex_lock(&root_mutex);
	if(!rootFile){
	
		// Open a root file to hold our output histogram(s)
		rootFile = new TFile(rootFileName,"RECREATE");
		if(!rootFile->IsOpen()){
			cerr<<"Unable to open ROOT output file!"<<endl;
			return UNKNOWN_ERROR;
		}

		// Create a histogram of dp over p vs. p
		dp_over_p_vs_p = new TH2D("dp_over_p_vs_p","#deltap/p vs. p", 100, 0.0, 7.0, 100, -0.2, 0.2);
	}

	pthread_mutex_unlock(&root_mutex);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_HDParSim::brun(JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_HDParSim::evnt(JEventLoop *loop, int eventnumber)
{
	// Fill histograms here

	// Smear momentum. On input, "mom" has the true values and on
	// output, it contains the smeared values.
	//res->Smear(PiPlus, mom);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_HDParSim::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_HDParSim::fini(void)
{
	// Close ROOT output file
	pthread_mutex_lock(&root_mutex);
	if(rootFile){
		rootFile->Write();
		delete rootFile;
		rootFile=NULL;
	}
	pthread_mutex_unlock(&root_mutex);

	return NOERROR;
}

