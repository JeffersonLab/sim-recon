// $Id$
//
//    File: DReaction_factory_p2pi_hists.cc
// Created: Wed Mar 11 20:34:14 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//


#include "DReaction_factory_p2pi_hists.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_p2pi_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("p2pi_pmiss"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p2pi_pmiss Reaction Steps ****************************************************/

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	// g, p -> pi+, pi- (p)
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(Gamma);
        locReactionStep->Set_TargetParticleID(Proton);
        locReactionStep->Add_FinalParticleID(PiPlus);
        locReactionStep->Add_FinalParticleID(PiMinus);
        locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p2pi_pmiss Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 2.004 ns)
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p2pi_pmiss Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// POCA cut on all tracks
        locReaction->Add_AnalysisAction(new DCutAction_AllVertexZ(locReaction, 50., 80.));

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	
	// Custom histograms for p2pi (no KinFit cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Require KinFit converges
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.0)); //require kinematic fit converges

	// Custom histograms for p2pi (KinFit converges)
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, false, "KinFitConverge_Measured"));

	// Require KinFit FOM > 0.1
        locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.1));

        // Custom histograms for p2pi (KinFit FOM > 10%)
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, false, "KinFitCut10_Measured"));
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, true, "KinFitCut10_KinFit"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory


	/**************************************************** p2pi_preco Reaction Steps ****************************************************/

	locReaction = new DReaction("p2pi_preco"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// g, p -> pi+, pi- ,p
        locReactionStep = new DReactionStep();
        locReactionStep->Set_InitialParticleID(Gamma);
        locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(Proton);
        locReactionStep->Add_FinalParticleID(PiPlus);
        locReactionStep->Add_FinalParticleID(PiMinus);
        locReaction->Add_ReactionStep(locReactionStep);
        dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p2pi_preco Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_NoFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	// Recommended: Enable ROOT TTree output for this DReaction
        locReaction->Enable_TTreeOutput("tree_p2pi_hists.root"); //string is file name (must end in ".root"!!): doen't need to be unique, feel free to change

	/**************************************************** p2pi_preco Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// POCA cut on all tracks
        locReaction->Add_AnalysisAction(new DCutAction_AllVertexZ(locReaction, 50., 80.));
	
	// Custom histograms for p2pi (no KinFit cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, false, "NoKinFit_Measured"));
	locReaction->Add_AnalysisAction(new DCustomAction_p2pi_unusedHists(locReaction, false, "NoKinFit_Measured"));

	// Kinematic Fit Results
	//locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Require KinFit converges
	//locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.0)); //require kinematic fit converges

	// Custom histograms for p2pi (KinFit converges)
        //locReaction->Add_AnalysisAction(new DCustomAction_p2pi_hists(locReaction, false, "KinFitConverge_Measured"));
	//locReaction->Add_AnalysisAction(new DCustomAction_p2pi_unusedHists(locReaction, false, "KinFitConverge_Measured"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p2pi_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

