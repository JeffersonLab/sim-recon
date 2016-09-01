// $Id$
//
//    File: DReaction_factory_p2k_hists.cc
// Created: Wed Mar 11 20:34:14 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DReaction_factory_p2k_hists.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_p2k_hists::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_p2k_hists::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("p2k_preco"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p2k_preco Reaction Steps ****************************************************/

	//locReaction = new DReaction("p2k_preco"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// g, p -> k+, k- ,p
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(Gamma);
        locReactionStep->Set_TargetParticleID(Proton);
        locReactionStep->Add_FinalParticleID(KPlus);
        locReactionStep->Add_FinalParticleID(KMinus);
        locReactionStep->Add_FinalParticleID(Proton);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p2k_preco Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-"); // boolean-AND of skims

	//locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)
	locReaction->Set_MaxExtraGoodTracks(4);

	/**************************************************** p2k_preco Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.75, KPlus, SYS_TOF)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, KPlus, SYS_BCAL)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, KPlus, SYS_FCAL)); //cut at delta-t +/- 1.0 //false: measured data

        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.75, KMinus, SYS_TOF)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, KMinus, SYS_BCAL)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, KMinus, SYS_FCAL)); //cut at delta-t +/- 1.0 //false: measured data

        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_TOF)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL)); //cut at delta-t +/- 1.0 //false: measured data
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL)); //cut at delta-t +/- 1.0 //false: measured data

	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// Custom histograms for p2k (no KinFit cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p2k_hists(locReaction, false, "NoKinFit_Measured"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

	// Kinematics fit
	// locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); // 5% confidence level cut on pull histograms only
	// locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); // -1.0 confidence level cut // require kinematic fit converges

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p2k_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

