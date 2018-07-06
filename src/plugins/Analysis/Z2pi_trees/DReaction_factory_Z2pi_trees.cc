// $Id$
//  DReaction_factory_Z2pi_trees, modeled after  DReaction_factory_p2pi_trees
//
//    File: DReaction_factory_p2pi_trees.cc
// Created: Wed Mar 29 16:34:58 EDT 2017
// Creator: elton (on Linux ifarm1401.jlab.org 3.10.0-327.el7.x86_64 x86_64)
// This DReation is geared toward gamma Z -> pi+ pi- Z, i.e. we do not expect to see the recoil nucleus.
//


// #include "DCustomAction_dEdxCut_Z2pi.h"
// #include "DCustomAction_Z2pi_trees.h"
// #include "DCustomAction_Z2pi_cuts.h"
// #include "DCustomAction_Z2pi_unusedHists.h"
#include "DReaction_factory_Z2pi_trees.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_Z2pi_trees::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_Z2pi_trees::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
  {
	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/************************************************** Z2pi_trees Reaction Definition *************************************************/

	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = NULL; //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	bool unused = false;
	locReaction = new DReaction("Z2pi_trees"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
	//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h
	// g, Z -> pi+, pi- ,Z   - assume lead for now
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Pb208);
	locReactionStep->Add_FinalParticleID(Pb208,true);   // recoil missing
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak


	/**************************************************** Z2pi_trees Control Settings ****************************************************/

	// Event Store
	// locReaction->Set_EventStoreSkims("q+,q-"); // boolean-AND of skims

	// Kinematic Fit
	//locReaction->Set_KinFitType(d_NoFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	// locReaction->Set_KinFitType(d_P4Fit); //simultaneously constrain apply four-momentum conservation, invariant masses, NO vertex

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	  // Outdated version for RF selection locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
	locReaction->Set_NumPlusMinusRFBunches(0);

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
		// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
	locReaction->Set_MaxExtraGoodTracks(4);

	// Highly Recommended: Enable ROOT TTree output for this DReaction
	// string is file name (must end in ".root"!!): doen't need to be unique, feel free to change
	locReaction->Enable_TTreeOutput("tree_Z2pi_trees.root", false); //true/false: do/don't save unused hypotheses

	/************************************************** Z2pi_trees Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass of proton
        locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, 20000, 60000));

	/**************************************************** Z2pi_trees Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
	//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination
	//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h
	
	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, false));
	/*locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));*/
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	// locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_Z2pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, false, "PostPIDCuts"));

	// Custom histograms for Z2pi  // omit for now since reaction is on heavy target
	// locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_hists(locReaction, false, "TimingCut_Measured"));

	//MISSING MASS
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, 30000, 50000, "PreKinFit"));
	// locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.005, "PreKinFit"));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, 30000,50000, "PreKinFit"));

	// Require KinFit converges
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.0)); //require kinematic fit converges

	if(unused) {
	  // Custom cuts (can be applied in TSelector)
	  // locReaction->Add_AnalysisAction(new DCutAction_ProtonPiPlusdEdx(locReaction, 2.2, true)); //select p/pi+ above/below 2.2, //true/false: cut all/no proton candidates above p = 1 GeV/c
	  locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, 30000,50000));
	  // locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_cuts(locReaction, false));

	  // Diagnostics for unused tracks and showers with final selection (only useful when analyzing EVIO data)
	  //comment locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_unusedHists(locReaction, false, "KinCut_Measured"));
	}


	// Custom histograms (after kinematic cuts)
	// locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_trees(locReaction, false, "KinCut_Measured"));

	// 2PI
	deque<Particle_t> loc2piPIDs;  loc2piPIDs.push_back(PiPlus);  loc2piPIDs.push_back(PiMinus);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, loc2piPIDs, false, 400, 0.2, 0.6, "2pi Mass"));

	if(unused)
	{
		// Custom cuts (can be applied in TSelector)
		// locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_cuts(locReaction, false));

		// Diagnostics for unused tracks and showers with final selection (only useful when analyzing EVIO data)
		//comment locReaction->Add_AnalysisAction(new DCustomAction_Z2pi_unusedHists(locReaction, false, "Unused"));
	}   

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_Z2pi_trees::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

