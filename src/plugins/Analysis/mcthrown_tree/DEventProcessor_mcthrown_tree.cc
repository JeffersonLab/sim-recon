// $Id$
//
// File: DEventProcessor_mcthrown_tree.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_mcthrown_tree.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"
{
	void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new DEventProcessor_mcthrown_tree());
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_mcthrown_tree::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_mcthrown_tree::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	const DEventWriterROOT* locEventWriterROOT = NULL;
	locEventLoop->GetSingle(locEventWriterROOT);
	locEventWriterROOT->Create_ThrownTree(locEventLoop, "tree_thrown.root");

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_mcthrown_tree::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	const DEventWriterROOT* locEventWriterROOT = NULL;
	locEventLoop->GetSingle(locEventWriterROOT);
	locEventWriterROOT->Fill_ThrownTree(locEventLoop);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_mcthrown_tree::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_mcthrown_tree::fini(void)
{
	return NOERROR;
}

