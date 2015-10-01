// $Id$
//
//    File: DReaction_factory_p3pi_hists.cc
// Created: Wed Mar 11 20:34:22 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//


#include "DReaction_factory_p3pi_hists.h"
#include "DCustomAction_HistOmegaVsMissProton.h"
#include "DCustomAction_dEdxCut.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_p3pi_hists::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction;

	double minPi0FCAL = 0.09;
	double maxPi0FCAL = 0.14;
	double minPi0BCAL = 0.11;
	double maxPi0BCAL = 0.16;
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
	locReactionStep->Set_ApplyKinFitMassConstraintOnInitialParticleFlag(false);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** p3pi_preco_2FCAL Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p3pi_preco_2FCAL Analysis Actions ****************************************************/

	// Require 2 photons in FCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 2));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Photon, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Photon, SYS_FCAL)); //false: measured data

	//Kinematics Pre-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Pre-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Pre-Pi0Cut"));

	// Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0FCAL, maxPi0FCAL));

	//Kinematics Post-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Post-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Post-Pi0Cut"));

	// Custom histograms for p3pi (after Pi0 mass cut)
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// dE/dx Cut (after custom action since it does dE/dx studies)
	locReaction->Add_AnalysisAction(new DCutAction_ProtonPiPlusdEdx(locReaction, 2.2, false)); //select p/pi+ above/below 2.2, //true/false: cut all/no proton candidates above p = 1 GeV/c

	//	Missing Pt
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingTransverseMomentum(locReaction, false, 500, 0.0, 1.0));

	//	Missing Mass Squared (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 600, -0.06, 0.06));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.005));

	// Omega Mass (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 500, 0.4, 1.4, "Omega_Kinfit"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, true, 0.7, 0.9));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_FCAL-BCAL Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_FCAL-BCAL"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco FCAL-BCAL Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p3pi_preco FCAL-BCAL Analysis Actions ****************************************************/

	// Require 1 photon in FCAL and 1 photon in BCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 1));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Photon, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Photon, SYS_FCAL)); //false: measured data

	//Kinematics Pre-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Pre-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Pre-Pi0Cut"));

	// Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0FCAL_BCAL, maxPi0FCAL_BCAL));

	//Kinematics Post-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Post-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Post-Pi0Cut"));

	// Custom histograms for p3pi (after Pi0 mass cut)
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// dE/dx Cut (after custom action since it does dE/dx studies)
	locReaction->Add_AnalysisAction(new DCutAction_ProtonPiPlusdEdx(locReaction, 2.2, false)); //select p/pi+ above/below 2.2, //true/false: cut all/no proton candidates above p = 1 GeV/c

	//	Missing Pt
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingTransverseMomentum(locReaction, false, 500, 0.0, 1.0));

	//	Missing Mass Squared (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 600, -0.06, 0.06));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.005));

	// Omega Mass (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 500, 0.4, 1.4, "Omega_Kinfit"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, true, 0.7, 0.9));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_2BCAL Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_2BCAL"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco 2BCAL Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4AndVertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/**************************************************** p3pi_preco 2BCAL Analysis Actions ****************************************************/

	// Require 1 photon in FCAL and 1 photon in BCAL
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 0));

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Photon, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Photon, SYS_FCAL)); //false: measured data

	//Kinematics Pre-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Pre-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Pre-Pi0Cut"));

	// Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0BCAL, maxPi0BCAL));

	//Kinematics Post-Pi0Cut
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction, "Post-Pi0Cut"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Post-Pi0Cut"));

	// Custom histograms for p3pi (after Pi0 mass cut)
	locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// dE/dx Cut (after custom action since it does dE/dx studies)
	locReaction->Add_AnalysisAction(new DCutAction_ProtonPiPlusdEdx(locReaction, 2.2, false)); //select p/pi+ above/below 2.2, //true/false: cut all/no proton candidates above p = 1 GeV/c

	//	Missing Pt
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingTransverseMomentum(locReaction, false, 500, 0.0, 1.0));

	//	Missing Mass Squared (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 600, -0.06, 0.06));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.005));

	// Omega Mass (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 500, 0.4, 1.4, "Omega_Kinfit"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, true, 0.7, 0.9));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_preco_any Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_preco_any"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"
	locReaction->Add_ReactionStep(dReactionStepPool[0]);
	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_preco_any Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_VertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	/*********************************************** p3pi_preco_any Pre-Combo Custom Cuts ***********************************************/

	// Loose Pi0 Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(Pi0, 0.0, 0.3);

	// Loose omega Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(omega, 0.5, 1.1);

	// Loose missing mass squared cut, applied just after creating the combination (before saving it)
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.06, 0.06));

	/**************************************************** p3pi_preco_any Analysis Actions ****************************************************/

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Photon, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Photon, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// Kinematic Fit: Vertex-Only Fit
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
//	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); //cut +/- 5 sigma
	locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7)); //cut +/- 5 sigma
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Pi0
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 600, 0.0, 0.3, "Pi0"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.0917604, 0.169925));

	//	Missing Mass Squared (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 600, -0.06, 0.06));
	locReaction->Add_AnalysisAction(new DCutAction_MissingMassSquared(locReaction, false, -0.01, 0.006));

	// Omega Mass (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, false, 600, 0.5, 1.1, "Omega_Kinfit"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, false, 0.71, 0.85));

	// Kinematics of final selection
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false, "Final")); //false: fill histograms with measured particle data

	_data.push_back(locReaction); //Register the DReaction with the factory



	/**************************************************** p3pi_pmiss_any Reaction Steps ****************************************************/

	locReaction = new DReaction("p3pi_pmiss_any"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// g, p -> omega, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	locReaction->Add_ReactionStep(dReactionStepPool[1]);
	locReaction->Add_ReactionStep(dReactionStepPool[2]);

	/**************************************************** p3pi_pmiss_any Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_VertexFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_MaxPhotonRFDeltaT(0.5*4.008); //beam bunches are every 4.008 ns, (2.004 should be minimum cut value)

	locReaction->Set_MaxExtraGoodTracks(1); //not ideal
	locReaction->Set_MaxNumBeamPhotonsInBunch(1); //not ideal: throws away a lot of signal

	/*********************************************** p3pi_pmiss_any Pre-Combo Custom Cuts ***********************************************/

	// Loose Pi0 Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(Pi0, 0.0, 0.3);

	// Loose omega Cut, Applied during Blueprint Construction
	locReaction->Set_InvariantMassCut(omega, 0.5, 1.1);

	// Loose missing mass squared cut, applied just after creating the combination (before saving it)
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.1, 2.56));

	/**************************************************** p3pi_pmiss_any Analysis Actions ****************************************************/

	// PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Photon, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Photon, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));

	// Kinematic Fit: Vertex-Only Fit
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
//	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 5.73303E-7)); //cut +/- 5 sigma
	locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7)); //cut +/- 5 sigma
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Pi0
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 600, 0.0, 0.3, "Pi0"));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.0917604, 0.169925));

	// Omega vs missing proton
	locReaction->Add_AnalysisAction(new DCustomAction_HistOmegaVsMissProton(locReaction));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, false, 0.65, 0.9));
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, false, 0.71, 0.85));

	//	Missing Mass Squared (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1064, -0.1, 2.56));

	// Omega Mass (Hist and Cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 600, 0.5, 1.1, "Omega_Kinfit"));

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

