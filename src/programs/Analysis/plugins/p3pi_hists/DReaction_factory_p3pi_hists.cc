// $Id$
//
//    File: DReaction_factory_p3pi_hists.cc
// Created: Wed Mar 11 20:34:22 EDT 2015
// Creator: jrsteven (on Linux halldw1.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//


#include "DReaction_factory_p3pi_hists.h"

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

        // g, p -> omega. p
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
        locReactionStep->Set_ApplyKinFitMassConstraintOnInitialParticleFlag(false);
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

	// POCA cut on all tracks
        locReaction->Add_AnalysisAction(new DCutAction_AllVertexZ(locReaction, 50., 80.));

        // Require 2 photons in FCAL
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 2));

        // Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
        locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0FCAL, maxPi0FCAL));

	// Custom histograms for p3pi (after Pi0 mass cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// Kinematics of final selection
        locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

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

	// POCA cut on all tracks
        locReaction->Add_AnalysisAction(new DCutAction_AllVertexZ(locReaction, 50., 80.));

        // Require 1 photon in FCAL and 1 photon in BCAL
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 1));

        // Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
        locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0FCAL_BCAL, maxPi0FCAL_BCAL));

	// Custom histograms for p3pi (after Pi0 mass cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// Kinematics of final selection
        locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

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

	// POCA cut on all tracks
        locReaction->Add_AnalysisAction(new DCutAction_AllVertexZ(locReaction, 50., 80.));

        // Require 1 photon in FCAL and 1 photon in BCAL
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_Pi0Cuts(locReaction, false, 0));

        // Custom histograms for p3pi (no KinFit cut)
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false,500,0.,1., "NoKinFit_Measured"));
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "NoKinFit_Measured"));

	// Pi0 mass cut
        locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, minPi0BCAL, maxPi0BCAL));

	// Custom histograms for p3pi (after Pi0 mass cut)
        locReaction->Add_AnalysisAction(new DCustomAction_p3pi_hists(locReaction, false, "CutPi0_Measured"));

	// Kinematics of final selection
        locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: fill histograms with measured particle data

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

