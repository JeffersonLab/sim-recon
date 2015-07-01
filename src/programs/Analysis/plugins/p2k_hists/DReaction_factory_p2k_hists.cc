// $Id$
//
//    File: DReaction_factory_p2k_hists.cc
// Created: Wed Mar 11 20:34:14 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DReaction_factory_p2k_hists.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_p2k_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("p2k_pmiss"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p2k_preco Reaction Steps ****************************************************/

	locReaction = new DReaction("p2k_preco"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

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

	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p2k_preco Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Unknown, SYS_TOF)); //false: measured data //Unknown: All PIDs
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 10.0, Unknown, SYS_BCAL)); //false: measured data //Unknown: All PIDs
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 10.0, Unknown, SYS_FCAL)); //false: measured data //Unknown: All PIDs

	// Custom histograms for p2k (no KinFit cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p2k_hists(locReaction, false, "NoKinFit_Measured"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

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

