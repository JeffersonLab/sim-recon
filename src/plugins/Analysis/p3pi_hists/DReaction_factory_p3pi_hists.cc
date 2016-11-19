// $Id$
//
//    File: DReaction_factory_p3pi_hists.cc
// Created: Wed Mar 11 20:34:22 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//


#include "DReaction_factory_p3pi_hists.h"
#include "DCustomAction_HistOmegaVsMissProton.h"
#include "DCustomAction_dEdxCut_p3pi.h"
#include "DCustomAction_CutExtraPi0.h"
#include "DCustomAction_CutExtraTrackPID.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_p3pi_hists::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// init
//------------------
jerror_t DReaction_factory_p3pi_hists::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction;

	double minPi0FCAL_BCAL = 0.11;
	double maxPi0FCAL_BCAL = 0.16;

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** p3pi_preco_2FCAL Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_2FCAL"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// g, p -> omega, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton); 
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	// omega -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	// pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p3pi_preco_2FCAL Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-,pi0"); // boolean-AND of skims

	// KINFIT
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod);

	/************************************************** p3pi_preco_2FCAL Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 0.1));

	/**************************************************** p3pi_preco_2FCAL Analysis Actions ****************************************************/

	// Require 2 photons in FCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 2));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_p3pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// MASSES
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit"));

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); // confidence level cut //+/- 5 sigma

	// MASSES, POST-KINFIT
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1, "PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit_PostKinFitCut"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_FCAL-BCAL Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_FCAL-BCAL"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco_FCAL-BCAL Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-,pi0"); // boolean-AND of skims

	// KINFIT
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod);

	/************************************************** p3pi_preco_FCAL-BCAL Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 0.1));

	/**************************************************** p3pi_preco FCAL-BCAL Analysis Actions ****************************************************/

	// Require 1 photon in FCAL and 1 photon in BCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 1));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_p3pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// MASSES
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit"));

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); // confidence level cut //+/- 5 sigma

	// MASSES, POST-KINFIT
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1, "PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit_PostKinFitCut"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_2BCAL Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_2BCAL"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco_2BCAL Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-,pi0"); // boolean-AND of skims

	// KINFIT
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod);

	/************************************************** p3pi_preco_2BCAL Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 0.1));

	/**************************************************** p3pi_preco_2BCAL Analysis Actions ****************************************************/

	// Require 1 photon in FCAL and 1 photon in BCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 0));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_p3pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// MASSES
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit"));

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); // confidence level cut //+/- 5 sigma

	// MASSES, POST-KINFIT
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1, "PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit_PostKinFitCut"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_any Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_any"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco_any Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-,pi0"); // boolean-AND of skims

	// KINFIT
	//locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod);

	/*********************************************** p3pi_preco_any Pre-Combo Custom Cuts ***********************************************/

	// Loose Pi0 Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);

	// Loose omega Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);

	// Loose missing mass squared cut, applied just after creating the combination (before saving it)
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 0.1));

	/**************************************************** p3pi_preco_any Analysis Actions ****************************************************/

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_p3pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// Pi0
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0FCAL_BCAL, maxPi0FCAL_BCAL));

	// Custom histograms for p3pi
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false));

	// MISSING MASS
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.005));

	// OMEGA
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_any_kinfit Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_any_kinfit"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco_any_kinfit Control Settings ****************************************************/

	// Event Store
	locReaction->Set_EventStoreSkims("2q+,q-,pi0"); // boolean-AND of skims

	// KINFIT
	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
	locReaction->Set_MaxPhotonRFDeltaT(0.5*dBeamBunchPeriod);

	/*********************************************** p3pi_preco_any_kinfit Pre-Combo Custom Cuts ***********************************************/

	// Loose Pi0 Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);

	// Loose omega Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(omega, 0.4, 1.2);

	// Loose missing mass squared cut, applied just after creating the combination (before saving it)
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 0.1));

	/**************************************************** p3pi_preco_any_kinfit Analysis Actions ****************************************************/

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.5, Proton, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Proton, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiPlus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_BCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, PiMinus, SYS_FCAL));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 3.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_p3pi(locReaction, false)); //false: focus on keeping signal
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// MASSES
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit"));

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); // confidence level cut //+/- 5 sigma

	// MASSES, POST-KINFIT
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 850, 0.05, 0.22, "Pi0_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1000, -0.1, 0.1, "PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega_PostKinFitCut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_KinFit_PostKinFitCut"));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_p3pi_hists::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}
