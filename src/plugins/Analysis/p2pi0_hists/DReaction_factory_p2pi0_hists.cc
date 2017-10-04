// $Id$
//
//    File: DReaction_factory_p2pi0_hists.cc
// Created: Tue May 19 07:56:15 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//


#include "DReaction_factory_p2pi0_hists.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_p2pi0_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("p2pi0_hists"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p2pi0_hists Reaction Steps ****************************************************/

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	// g, p -> pi0, pi0, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	// pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Set_KinFitConstrainInitMassFlag(false);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	// pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Set_KinFitConstrainInitMassFlag(false);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p2pi0_hists Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 2.004 ns)
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p2pi0_hists Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, false));
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Unknown, SYS_TOF)); //false: measured data //Unknown: All PIDs
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 10.0, Unknown, SYS_BCAL)); //false: measured data //Unknown: All PIDs
        locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 10.0, Unknown, SYS_FCAL)); //false: measured data //Unknown: All PIDs

        // Custom histograms for p2pi0 (no KinFit cut)
        locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi0_hists(locReaction, false, "NoKinFit_Measured"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p2pi0_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

