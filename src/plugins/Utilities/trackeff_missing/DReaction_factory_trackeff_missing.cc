// $Id$
//
//    File: DReaction_factory_trackeff_missing.cc
// Created: Wed Feb 25 08:58:19 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#include "DReaction_factory_trackeff_missing.h"
#include "DCustomAction_TrackingEfficiency.h"
#include "DCustomAction_CutExtraPi0.h"
#include "DCustomAction_CutExtraShowers.h"
#include "DCustomAction_dEdxCut_trackeff.h"

//------------------
// brun
//------------------
jerror_t DReaction_factory_trackeff_missing::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_trackeff_missing::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// Make as many DReaction objects as desired
	DReactionStep* locReactionStep = NULL;
	DReaction* locReaction = new DReaction("TrackEff_MissingProton"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	/**************************************************** TrackEff_MissingProton Reaction Steps ****************************************************/

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi+, pi-, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingProton Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
	//locReaction->Set_KinFitType(d_P4Fit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_KinFitType(d_P4Fit);

	// Highly Recommended: When generating particle combinations, reject all photon candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinPhotonPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all charged track candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinChargedPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
		// Current (09/26/2014): "Good" tracks have a detector-hit match, and tracking FOM > 0.0027 (+/- 3 sigma). 
		// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingProton Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, 0.5, 1.4));

	/**************************************************** TrackEff_MissingProton Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));		

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Kinematic Fit: Vertex-Only Fit
	//locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7, 44.0, 85.0)); //cut +/- 5 sigma, vertex-z between 44 & 85

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	//_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** TrackEff_MissingPiMinus Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiMinus"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi+, p, (pi-)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus, true); //true: pi- missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiMinus Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
	//locReaction->Set_KinFitType(d_P4Fit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_KinFitType(d_P4Fit);

	// Highly Recommended: When generating particle combinations, reject all photon candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinPhotonPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all charged track candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinChargedPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
		// Current (09/26/2014): "Good" tracks have a detector-hit match, and tracking FOM > 0.0027 (+/- 3 sigma).
		// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiMinus Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiMinus Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

//	locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction, false, 3.0, 6.0));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Kinematic Fit: Vertex-Only Fit
	//locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7, 44.0, 85.0)); //cut +/- 5 sigma, vertex-z between 44 & 85

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6, "PostPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	_data.push_back(locReaction); //Register the DReaction with the factory


	/**************************************************** TrackEff_MissingPiPlus Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiPlus"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi-, p, (pi+)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus, true); //true: pi+ missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiPlus Control Settings ****************************************************/

	// Recommended: Type of kinematic fit to perform (default is d_NoFit)
		//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DKinFitResults.h
	//locReaction->Set_KinFitType(d_P4Fit); //simultaneously constrain apply four-momentum conservation, invariant masses, and common-vertex constraints
	locReaction->Set_KinFitType(d_P4Fit);

	// Highly Recommended: When generating particle combinations, reject all photon candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinPhotonPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all charged track candidates with a PID confidence level < 5.73303E-7 (+/- 5-sigma)
	// locReaction->Set_MinChargedPIDFOM(5.73303E-7);

	// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch (delta_t > 1.002 ns)
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value

	// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
		// Current (09/26/2014): "Good" tracks have a detector-hit match, and tracking FOM > 0.0027 (+/- 3 sigma).
		// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiPlus Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiPlus Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

//	locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction, false, 3.0, 6.0));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Kinematic Fit: Vertex-Only Fit
	//locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7, 44.0, 85.0)); //cut +/- 5 sigma, vertex-z between 44 & 85

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6, "PostPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/*************************************************************************************************************/
	/**************************************************** 3pi ****************************************************/
	/*************************************************************************************************************/

	/**************************************************** TrackEff_MissingProton_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	deque<Particle_t> locToIncludePIDs;
	locToIncludePIDs.push_back(PiPlus); locToIncludePIDs.push_back(PiMinus); locToIncludePIDs.push_back(Pi0); 

	//Example: g, p -> pi+, pi-, pi0 , (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingProton_3pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingProton_3pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, 0.5, 1.4));

	/**************************************************** TrackEff_MissingProton_3pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

//	locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction, false, 3.0, 6.0));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));
	
	// Missing Mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 600, 0.5, 1.7, "3pi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));

	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 600, 0.5, 1.7, "3pi_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostExtraPi0"));
	
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.115, 0.155));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 600, 0.5, 1.7, "3pi_PostPi0Mass"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostPi0Mass"));
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true, "PostPi0Mass")); //true: fill histograms with kinematic-fit particle data
	
	//locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 1.2, 1.7));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "Post3piMass"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	//_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** TrackEff_MissingPiPlus_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiPlus_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Example: g, p -> (pi+), pi-, pi0 , p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus, true); // true: pi+ missing
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Proton); 
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiPlus_3pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiPlus_3pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiPlus_3pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));
	
	// Missing Mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));

	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4, "PostExtraPi0"));
	
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.115, 0.155));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi_PostPi0Mass"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4, "Post3piMass"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	//_data.push_back(locReaction); //Register the DReaction with the factory

	/**************************************************** TrackEff_MissingPiMinus_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiMinus_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Example: g, p -> pi+, (pi-), pi0 , p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus, true); // true: pi- missing
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiMinus_3pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiMinus_3pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	locReaction->Set_InvariantMassCut(Pi0, 0.05, 0.22);

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiMinus_3pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));
	
	// Missing Mass
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));

	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4, "PostExtraPi0"));
	
	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.115, 0.155));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, 0, locToIncludePIDs, true, 300, 0.5, 1.1, "3pi_PostPi0Mass"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 500, -0.3, 0.4, "Post3piMass"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	//_data.push_back(locReaction); //Register the DReaction with the factory

	/*************************************************************************************************************/
	/**************************************************** 4pi ****************************************************/
	/*************************************************************************************************************/

	/**************************************************** TrackEff_MissingProton_4pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton_4pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> 2pi+, 2pi-, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingProton_4pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingProton_4pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, 0.5, 1.4));

	/**************************************************** TrackEff_MissingProton_4pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

//	locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction, false, 3.0, 6.0));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Kinematic Fit: Vertex-Only Fit
	//locReaction->Add_AnalysisAction(new DCutAction_OneVertexKinFit(locReaction, 5.73303E-7, 44.0, 85.0)); //cut +/- 5 sigma, vertex-z between 44 & 85

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));

	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostExtraPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	_data.push_back(locReaction); //Register the DReaction with the factory


	/**************************************************** TrackEff_MissingPiPlus_4pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiPlus_4pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> pi+, (pi+), 2pi-, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus, true); //true: piplus missing
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton); 
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiPlus_4pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiPlus_4pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiPlus_4pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6, "PostExtraPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	_data.push_back(locReaction); //Register the DReaction with the factory


	/**************************************************** TrackEff_MissingPiMinus_4pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiMinus_4pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Required: DReactionSteps to specify the channel and decay chain you want to study
		//Particles are of type Particle_t, an enum defined in sim-recon/src/libraries/include/particleType.h

	//Example: g, p -> 2pi+, pi-, (pi-) p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus); 
	locReactionStep->Add_FinalParticleID(PiMinus, true); //true: piminus missing
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton); 
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingPiMinus_4pi Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingPiMinus_4pi Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, -0.5, 0.6));

	/**************************************************** TrackEff_MissingPiMinus_4pi Analysis Actions ****************************************************/

	// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
		//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
		//Pre-defined actions can be found in ANALYSIS/DHistogramActions.h and ANALYSIS/DCutActions.h

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, Proton, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.8, Proton, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Proton, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiPlus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, PiMinus, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Track DOCA, etc.
	locReaction->Add_AnalysisAction(new DHistogramAction_TrackVertexComparison(locReaction));

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6));
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
	locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, 1100, -0.5, 0.6, "PostExtraPi0"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	_data.push_back(locReaction); //Register the DReaction with the factory

	/*************************************************************************************************************/
	/******************************************* omega -> pi0 gamma **********************************************/
	/*************************************************************************************************************/

	/**************************************************** TrackEff_MissingProton_Omega Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton_Omega"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//Example: g, p -> omega, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//omega -> pi0, gamma
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	/**************************************************** TrackEff_MissingProton_Omega Control Settings ****************************************************/

	locReaction->Set_KinFitType(d_P4Fit);
	locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); //should be minimum cut value
//	locReaction->Set_MaxExtraGoodTracks(1);

	/************************************************** TrackEff_MissingProton_Omega Pre-Combo Custom Cuts *************************************************/

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	// Example: Missing mass squared of proton
	locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, 0.5, 1.4));

	/**************************************************** TrackEff_MissingProton_Omega Analysis Actions ****************************************************/

	// Reject events with large unused shower energy
	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

//	locReaction->Add_AnalysisAction(new DCutAction_BeamEnergy(locReaction, false, 3.0, 6.0));		

	// Hist PID
	locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

	// PID Cuts
	locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Gamma, SYS_BCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_FCAL)); //false: measured data
	locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
	locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

	// Kinematic Fit Results
	locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only

	// Missing Mass Squared
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 300, 0.5, 1.1, "Omega"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4));

	locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));	
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, Pi0, false, 400, 0.05, 0.22, "Pi0_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 300, 0.5, 1.1, "Omega_PostExtraPi0"));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostExtraPi0"));

	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, Pi0, false, 0.115, 0.155));
	locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, omega, true, 300, 0.5, 1.1, "Omega_PostPi0Mass"));

	locReaction->Add_AnalysisAction(new DCutAction_InvariantMass(locReaction, omega, true, 0.757, 0.807));
	locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, 900, 0.5, 1.4, "PostOmegaMass"));

	// Kinematics
	locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

	// Tracking Efficiency
	locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

	//_data.push_back(locReaction); //Register the DReaction with the factory

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_trackeff_missing::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

