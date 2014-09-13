// $Id$
//
// File: DEventProcessor_monitoring_hists.cc
// Created: Thu Sep 28 11:38:03 EDT 2011
// Creator: pmatt (on Darwin swire-b241.jlab.org 8.4.0 powerpc)
//

#include "DEventProcessor_monitoring_hists.h"

// The executable should define the ROOTfile global variable. It will
// be automatically linked when dlopen is called.
extern TFile *ROOTfile;

// Routine used to create our DEventProcessor
extern "C"
{
	void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new DEventProcessor_monitoring_hists());
		app->AddFactoryGenerator(new DFactoryGenerator_monitoring_hists());
	}
} // "C"

//------------------
// init
//------------------
jerror_t DEventProcessor_monitoring_hists::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_monitoring_hists::brun(JEventLoop *locEventLoop, int runnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//Initialize Actions
	dHistogramAction_TrackMultiplicity.Initialize(locEventLoop);
	dHistogramAction_DetectedParticleKinematics.Initialize(locEventLoop);
	dHistogramAction_NumReconstructedObjects.Initialize(locEventLoop);
	dHistogramAction_DetectorStudies.Initialize(locEventLoop);

	if(!locMCThrowns.empty())
	{
		dHistogramAction_ThrownParticleKinematics.Initialize(locEventLoop);
		dHistogramAction_GenReconTrackComparison.Initialize(locEventLoop);
		dHistogramAction_ReconnedThrownKinematics.Initialize(locEventLoop);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_monitoring_hists::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//Fill reaction-independent histograms.
	dHistogramAction_TrackMultiplicity(locEventLoop);
	dHistogramAction_DetectedParticleKinematics(locEventLoop);
	dHistogramAction_NumReconstructedObjects(locEventLoop);
	dHistogramAction_DetectorStudies(locEventLoop);

	if(!locMCThrowns.empty())
	{
		dHistogramAction_ThrownParticleKinematics(locEventLoop);
		dHistogramAction_GenReconTrackComparison(locEventLoop);
		dHistogramAction_ReconnedThrownKinematics(locEventLoop);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_monitoring_hists::erun(void)
{
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventProcessor_monitoring_hists::fini(void)
{
	return NOERROR;
}

