// $Id$
//
//    File: DReaction_factory_pi0calib.cc
// Created: Tue Apr 28 21:19:40 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//


#include "DReaction_factory_pi0calib.h"
#include "DCustomAction_CutPhotonKin.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_pi0calib::init(void)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("excl_pi0calib"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** pi0calib Reaction Steps ****************************************************/

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//Example: pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Set_KinFitConstrainInitMassFlag(false);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** pi0calib Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
//	locReaction->Set_KinFitType(d_NoFit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_KinFitType(d_P4AndVertexFit);

	// Recommended: Enable ROOT TTree output for this DReaction
//        locReaction->Enable_TTreeOutput("tree_pi0calib.root"); //string is file name (must end in ".root"!!): doen't need to be unique, feel free to change

	/**************************************************** pi0calib Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
	//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
	//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	locReaction->Set_NumPlusMinusRFBunches(0);   // no accidentals!
	locReaction->Set_MaxExtraGoodTracks(1);
	//locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);   // Use default AnLib cut of M(gg) = 80 - 190 MeV
	
	// Require BCAL photons
	locReaction->Add_AnalysisAction(new DCustomAction_CutPhotonKin(locReaction));
	
	// Make some back-to-backness and other exclusivity cuts 
	locReaction->Add_AnalysisAction(new DCustomAction_p2gamma_cuts(locReaction, false));
	
	// Require kin fit CL > 1%
	locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, 0.01));

	_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_pi0calib::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

