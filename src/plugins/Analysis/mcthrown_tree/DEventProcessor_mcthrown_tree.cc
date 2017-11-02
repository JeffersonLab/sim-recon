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
// evnt
//------------------
jerror_t DEventProcessor_mcthrown_tree::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
	const DEventWriterROOT* locEventWriterROOT = NULL;
	locEventLoop->GetSingle(locEventWriterROOT);

	//This looks bad.  Really bad.  But relax, it's fine. 
	//This was previously in brun(), but brun() is ONLY CALLED BY ONE THREAD (no matter how many threads you run with!)
	//So this meant that the event writer wasn't set up for other threads!!
	//So, we call this here, every event, where every thread can see it. 
	//However, only the FIRST thread will actually create the tree
	//And, only the FIRST CALL for each thread will setup the event writer
	//Subsequent calls will auto-detect that everything is already done and bail early. 
	//Ugly, yes. But it works, and it's too late now to have each thread call brun().
	locEventWriterROOT->Create_ThrownTree(locEventLoop, "tree_thrown.root");

	// only keep generated events which hit a tagger counter
	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "TAGGEDMCGEN");

	// skip events where generated beam photon did not hit TAGM or TAGH counter
	if(dTagCheck && locBeamPhotons.empty())
		return NOERROR;

	locEventWriterROOT->Fill_ThrownTree(locEventLoop);

	return NOERROR;
}

