// $Id$
//
// File: DEventProcessor_b1pi_hists.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_b1pi_hists.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"
{
	void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new DEventProcessor_b1pi_hists());
		app->AddFactoryGenerator(new DFactoryGenerator_DReaction());
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_b1pi_hists::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_b1pi_hists::brun(JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_b1pi_hists::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_b1pi_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_b1pi_hists::fini(void)
{
	return NOERROR;
}

