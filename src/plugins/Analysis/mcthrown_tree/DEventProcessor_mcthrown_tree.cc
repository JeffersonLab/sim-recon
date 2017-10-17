// $Id$
//
// File: DEventProcessor_mcthrown_tree.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_mcthrown_tree.h"

#include <TAGGER/DTAGHGeometry.h>

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
	// require tagger hit for MCGEN beam photon by default to write event to TTree
	dTagCheck = true;
	gPARMS->SetDefaultParameter("MCTHROWN:TAGCHECK", dTagCheck);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_mcthrown_tree::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
	dNC_TAGH = DTAGHGeometry::kCounterCount;

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
	// only keep generated events which hit a tagger counter
        const DBeamPhoton* locBeamPhoton = NULL;
        locEventLoop->GetSingle(locBeamPhoton, "MCGEN");

	// skip events where generated beam photon did not hit TAGM or TAGH counter (ie. dCounter > 274)
	if(dTagCheck && locBeamPhoton->dCounter > dNC_TAGH)
		return NOERROR;

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

