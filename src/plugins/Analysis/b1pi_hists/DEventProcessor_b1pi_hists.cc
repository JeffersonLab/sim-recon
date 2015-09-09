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
		app->AddFactoryGenerator(new DFactoryGenerator_b1pi_hists());
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
	//Create Trees
	const DEventWriterROOT* locEventWriterROOT = NULL;
	locEventLoop->GetSingle(locEventWriterROOT);
	locEventWriterROOT->Create_DataTrees(locEventLoop);
	locEventWriterROOT->Create_ThrownTree("tree_b1pi_thrownmc.root");

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_b1pi_hists::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	//Triggers the analysis (is also automatically called by DEventWriterROOT::Fill_Trees())
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	//Output TTree
	const DEventWriterROOT* locEventWriterROOT = NULL;
	locEventLoop->GetSingle(locEventWriterROOT);
	locEventWriterROOT->Fill_DataTrees(locEventLoop, "b1pi_hists");
	locEventWriterROOT->Fill_ThrownTree(locEventLoop);

	//Do Miscellaneous Cuts
	bool locSaveEventFlag = false;
	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		const DAnalysisResults* locAnalysisResults = locAnalysisResultsVector[loc_i];
		if(locAnalysisResults->Get_Reaction()->Get_ReactionName() != "b1pi")
			continue;
		if(locAnalysisResults->Get_NumPassedParticleCombos() == 0)
			continue;
		locSaveEventFlag = true;
		break;
	}

	//Output REST File
	if(locSaveEventFlag)
	{
		vector<const DEventWriterREST*> locEventWriterRESTVector;
		locEventLoop->Get(locEventWriterRESTVector);
		locEventWriterRESTVector[0]->Write_RESTEvent(locEventLoop, "b1pi");
	}

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

