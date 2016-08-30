// $Id$
//
//    File: DReaction_factory_p4pi_hists.cc
// Created: Mon Aug 29 16:20:57 EDT 2016
// Creator: aaustreg (on Linux halld01.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//


#include "DReaction_factory_p4pi_hists.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_p4pi_hists::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_p4pi_hists::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = NULL; //create with a unique name for each DReaction object. CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/************************************************** p4pi_hists Reaction Definition *************************************************/

	locReaction = new DReaction("p4pi_hists");

	
	//Required: DReactionSteps to specify the channel and decay chain you want to study
	//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h
	
	// g, p -> pi+, pi-, pi+, pi-, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak


	/**************************************************** p4pi_hists Control Settings ****************************************************/

	// Highly Recommended: Set EventStore skim query (use with "eventstore" source)
	// This will skip creating particle combos for events that aren't in the skims you list
	// Query should be comma-separated list of skims to boolean-AND together
	locReaction->Set_EventStoreSkims("3q+,2q-"); //boolean-AND of skims

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
	//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DReaction.h
	//Options: d_NoFit (default), d_P4Fit, d_VertexFit, d_P4AndVertexFit
	//P4 fits automatically constrain decaying particle masses, unless they are manually disabled
	locReaction->Set_KinFitType(d_P4Fit);

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod); //should be minimum cut value


	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
	// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
	locReaction->Set_MaxExtraGoodTracks(4);

	// Highly Recommended: Enable ROOT TTree output for this DReaction
	// string is file name (must end in ".root"!!): doen't need to be unique, feel free to change
	// locReaction->Enable_TTreeOutput("tree_p4pi_hists.root", false); //true/false: do/don't save unused hypotheses

	/************************************************** p4pi_hists Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, -0.2, 0.2));

	/**************************************************** p4pi_hists Analysis Actions ****************************************************/

	
	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
	//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
	//Pre-defined actions can be found in ANALYSIS/DHistogramActions_*.h and ANALYSIS/DCutActions.h
	//If a histogram action is repeated, it should be created with a unique name (string) to distinguish them
	
	// HISTOGRAM PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// CUT PID
	// SYS_TOF, SYS_BCAL, SYS_FCAL, ...: DetectorSystem_t: Defined in libraries/include/GlueX.h
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL)); //false: measured data
	
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL)); //false: measured data
	
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL)); //false: measured data
	
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));


	// HISTOGRAM MASSES //false/true: measured/kinfit data
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, 0.04, 0.04, "PreKinFit"));
	
	std::deque<Particle_t> Four, Two, Zp, Zm;
	Four.push_back(PiPlus); Four.push_back(PiMinus); Four.push_back(PiPlus); Four.push_back(PiMinus);
	Two.push_back(PiPlus); Two.push_back(PiMinus);
	Zp.push_back(PiPlus); Zp.push_back(Proton);
	Zm.push_back(PiMinus); Zm.push_back(Proton);
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, Four, false, 500, 0.5, 2.5, "FourPi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, Two, false, 500, 0.0, 2.0, "TwoPi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_2DInvariantMass(locReaction, 0, Two, Two, false, 400, 0.2, 2.2, 400, 0.2, 2.2, "TwoPi_vs_TwoPi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_Dalitz(locReaction, 0, Two, Two, false, 600, 0, 3, 600, 0, 3, "Dalitz"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, Zp, false, 500, 1.0, 3.0, "ProtonPip"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, Zm, false, 500, 1.0, 3.0, "ProtonPim"));
	
	// KINEMATIC FIT
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); // -1.0 confidence level cut //require kinematic fit converges
	
	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "KinFit")); //true: kinematic-fit data
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));
	
	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p4pi_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

