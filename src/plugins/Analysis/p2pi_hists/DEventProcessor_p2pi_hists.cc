// $Id$
//
//    File: DEventProcessor_p2pi_hists.cc
// Created: Wed Mar 11 20:34:14 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DEventProcessor_p2pi_hists.h"

// Routine used to create our DEventProcessor

extern "C"
{
	void InitPlugin(JApplication *locApplication)
	{
		InitJANAPlugin(locApplication);
		locApplication->AddProcessor(new DEventProcessor_p2pi_hists()); //register this plugin
		locApplication->AddFactoryGenerator(new DFactoryGenerator_p2pi_hists()); //register the factory generator
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_p2pi_hists::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... create historgrams or trees ...
	// japp->RootUnLock();
	//

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_p2pi_hists::brun(jana::JEventLoop* locEventLoop, int locRunNumber)
{
	// This is called whenever the run number changes

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_p2pi_hists::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// locEventLoop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	//
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// locEventLoop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software

	//Recommended: Write surviving particle combinations (if any) to output ROOT TTree
                //If no cuts are performed by the analysis actions added to a DReaction, then this saves all of its particle combinations. 
                //The event writer gets the DAnalysisResults objects from JANA, performing the analysis. 
        // string is DReaction factory tag: will fill trees for all DReactions that are defined in the specified factory
        const DEventWriterROOT* locEventWriterROOT = NULL;
        locEventLoop->GetSingle(locEventWriterROOT);
        //locEventWriterROOT->Fill_DataTrees(locEventLoop, "p2pi_hists");

	// Get the analysis results for all DReactions. 
		//Getting these objects triggers the analysis, if it wasn't performed already. 
		//These objects contain the DParticleCombo objects that survived the DAnalysisAction cuts that were added to the DReactions
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_p2pi_hists::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_p2pi_hists::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

