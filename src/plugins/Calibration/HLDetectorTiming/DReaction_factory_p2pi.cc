// $Id$
//
//    File: DReaction_factory_p2pi.cc
// Created: Thu May  7 16:22:04 EDT 2015
// Creator: mstaib (on Linux gluon109.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//


#include "DReaction_factory_p2pi.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_p2pi::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("p2pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p2pi Reaction Steps ****************************************************/
 
	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi+, pi-, pi0, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton, false); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//optional: in a p4 kinematic fit, this will disable the constraint on the mass of the pi0
		//default is enabled, but you don't need to worry about setting false for the beam particle step
		//note that the omega and phi will always be false, regardless of what you set below
	//locReactionStep->Set_ApplyKinFitMassConstraintOnInitialParticleFlag(false);

	/**************************************************** p2pi Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
    locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all photon candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	//locReaction->Set_MinPhotonPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all charged track candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	//locReaction->Set_MinChargedPIDFOM(5.73303E-7);

	// Highly Recommended: When recording thrown/reconstructed matches in ROOT TTree, only register matches with match FOM > 5.73303E-7 (+/- 5-sigma)
	//locReaction->Set_MinThrownMatchFOMForROOT(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	//locReaction->Set_MaxPhotonRFDeltaT(0.5*2.004); //beam bunches are every 2.004 ns, (1.002 should be minimum cut value)

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
		// Current (09/26/2014): "Good" tracks have a detector-hit match, and tracking FOM > 0.0027 (+/- 3 sigma). 
		// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
	locReaction->Set_MaxExtraGoodTracks(0);

	// Recommended: Enable ROOT TTree output for this DReaction
	// locReaction->Enable_TTreeOutput("tree_p2pi.root"); //string is file name (must end in ".root"!!): doen't need to be unique, feel free to change

	/************************************************** p2pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	// Example: pi0 -> g, g cut
	// locReaction->Set_InvariantMassCut(Pi0, 0.08, 0.19);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	// locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 2.56));

	/**************************************************** p2pi Analysis Actions ****************************************************/

	
	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// PID & Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	//locReaction->Add_AnalysisAction(new DHistogramAction_TruePID(locReaction)); //momentum distributions of tracks with true/false PID (if thrown data available)

	// Transverse Momentum //Recommended for no-missing-particle reactions only!
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingTransverseMomentum(locReaction, false, 1000, 0.0, 1.0)); //false: fill histograms with measured particle data
    locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
    locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.03, 0.03));
	locReaction->Add_AnalysisAction(new DCutAction_TransverseMomentum(locReaction, 0.4)); //Max Missing Pt of 0.4 GeV

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	// locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.0)); //0% confidence level cut //require kinematic fit converges

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data
	// locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "KinFit")); //true: fill histograms with kinematic-fit particle data //"KinFit": unique name since action type is repeated
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));
	

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p2pi::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

