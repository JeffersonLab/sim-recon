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
#include "DCustomAction_CutNoDetectorHit.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_trackeff_missing::init(void)
{
	Define_LooseCuts();
	return NOERROR;
}

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

	// DEFINE CHANNELS:
	deque<DReaction*> locReactions;


	/**************************************************** TrackEff_MissingProton Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> pi+, pi-, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingPiMinus Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiMinus"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> pi+, p, (pi-)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus, true); //true: pi- missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingPiPlus Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiPlus"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> pi-, p, (pi+)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Proton);
	locReactionStep->Add_FinalParticleID(PiPlus, true); //true: pi+ missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingProton_4pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton_4pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

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
	locReactions.push_back(locReaction);

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
	locReactions.push_back(locReaction);

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
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingProton_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingProton_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> omega, (p)
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton, true); //true: proton missing
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//omega -> pi+, pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingPiMinus_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiMinus_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> omega, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//omega -> pi+, (pi-), pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus);
	locReactionStep->Add_FinalParticleID(PiMinus, true); //true: pi- missing
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);

	/**************************************************** TrackEff_MissingPiPlus_3pi Reaction Steps ****************************************************/

	locReaction = new DReaction("TrackEff_MissingPiPlus_3pi"); //needs to be a unique name for each DReaction object, CANNOT (!) be "Thrown"

	//g, p -> omega, p
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);
	locReactionStep->Add_FinalParticleID(omega);
	locReactionStep->Add_FinalParticleID(Proton);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//omega -> (pi+), pi-, pi0
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(omega);
	locReactionStep->Add_FinalParticleID(PiPlus, true); //true: pi+ missing
	locReactionStep->Add_FinalParticleID(PiMinus);
	locReactionStep->Add_FinalParticleID(Pi0);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak

	//pi0 -> g, g
	locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Pi0);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReactionStep->Add_FinalParticleID(Gamma);
	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep); //register so will be deleted later: prevent memory leak
	locReactions.push_back(locReaction);





	for(auto& locReaction : locReactions)
	{
		/**************************************************** Control Settings ****************************************************/

		// Recommended: Type of kinematic fit to perform (default is d_NoFit)
			//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DReaction.h
			//Options: d_NoFit (default), d_P4Fit, d_VertexFit, d_P4AndVertexFit
			//P4 fits automatically constrain decaying particle masses, unless they are manually disabled
		locReaction->Set_KinFitType(d_P4Fit); //d_P4AndVertexFit //No vertex: can't cut on kinfit conlev anyway, but could distort if alignment is bad
		locReaction->Set_KinFitUpdateCovarianceMatricesFlag(true);

		// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
		locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); // +/- 1 bunch for sideband subtraction

		/************************************************** Pre-Combo Custom Cuts *************************************************/

		Add_PreComboCuts(locReaction);

		// PID
		Add_PIDActions(locReaction);

		/**************************************************** Analysis Actions ****************************************************/

		// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
			//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination
			//Pre-defined actions can be found in ANALYSIS/DHistogramActions_*.h and ANALYSIS/DCutActions.h
			//If a histogram action is repeated, it should be created with a unique name (string) to distinguish them

		//FURTHER PID
		locReaction->Add_AnalysisAction(new DCustomAction_dEdxCut_trackeff(locReaction, true)); //true: focus on rejecting background
		locReaction->Add_AnalysisAction(new DCutAction_TrackFCALShowerEOverP(locReaction, false, 0.5)); //false: measured data //value: cut e+/e- below this, tracks above this
		locReaction->Add_AnalysisAction(new DCutAction_EachPIDFOM(locReaction, -9.9E9, true)); //cut particles with PID FOM = 0

		// SHOWER BACKGROUND
		locReaction->Add_AnalysisAction(new DCustomAction_CutExtraPi0(locReaction, 0.0775209, 0.188047));
		locReaction->Add_AnalysisAction(new DCustomAction_CutExtraShowers(locReaction, 0.5));

		//TRACK PURITY
		locReaction->Add_AnalysisAction(new DCutAction_MinTrackHits(locReaction, 10));

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, false, "PreKinFit");

		// KINEMATIC FIT
		locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
		locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); //require kinematic fit converges

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		locReaction->Add_AnalysisAction(new DCustomAction_CutNoDetectorHit(locReaction));
		Add_MassHistograms(locReaction, false, "PostKinFit");
		Add_MassHistograms(locReaction, true, "PostKinFit_KinFit");

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: fill histograms with kinematic-fit particle data

		// Tracking Efficiency
		locReaction->Add_AnalysisAction(new DCustomAction_TrackingEfficiency(locReaction, true));

		_data.push_back(locReaction); //Register the DReaction with the factory
	}

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

/************************************************************** ACTIONS AND CUTS **************************************************************/

void DReaction_factory_trackeff_missing::Define_LooseCuts(void)
{
	// Missing Mass Cuts
	dMissingMassCuts[Unknown] = pair<double, double>(-0.1, 0.1);
	dMissingMassCuts[Proton] = pair<double, double>(-0.5, 3.0);
	dMissingMassCuts[Neutron] = dMissingMassCuts[Proton];
        dMissingMassCuts[PiPlus] = pair<double, double>(-1.0, 1.0);
	dMissingMassCuts[PiMinus] = dMissingMassCuts[PiPlus];
	dMissingMassCuts[omega] = pair<double, double>(0.7, 0.9);

	// Invariant Mass Cuts: Mesons
	dInvariantMassCuts[Pi0] = pair<double, double>(0.08, 0.19);
	dInvariantMassCuts[KShort] = pair<double, double>(0.3, 0.7);
	dInvariantMassCuts[Eta] = pair<double, double>(0.3, 0.8);
	dInvariantMassCuts[omega] = pair<double, double>(0.74, 0.83);
	dInvariantMassCuts[EtaPrime] = pair<double, double>(0.6, 1.3);
	dInvariantMassCuts[phiMeson] = pair<double, double>(0.8, 1.2);
	dInvariantMassCuts[Jpsi] = pair<double, double>(1.5, 4.0);

	// Invariant Mass Cuts: Baryons
	dInvariantMassCuts[Lambda] = pair<double, double>(1.0, 1.2);
	dInvariantMassCuts[Sigma0] = pair<double, double>(1.1, 1.3);
	dInvariantMassCuts[SigmaPlus] = pair<double, double>(1.1, 1.3);
	dInvariantMassCuts[XiMinus] = pair<double, double>(1.1, 1.5);
	dInvariantMassCuts[Xi0] = pair<double, double>(1.1, 1.5);

	// Timing Cuts: Photon
	dPIDTimingCuts[Gamma][SYS_BCAL] = 0.4;
	dPIDTimingCuts[Gamma][SYS_FCAL] = 1.0;

	// Timing Cuts: Leptons
	dPIDTimingCuts[Electron][SYS_BCAL] = 0.6;
	dPIDTimingCuts[Electron][SYS_FCAL] = 1.4;
	dPIDTimingCuts[Electron][SYS_TOF] = 0.3;
	dPIDTimingCuts[Positron] = dPIDTimingCuts[Electron];
	dPIDTimingCuts[MuonMinus] = dPIDTimingCuts[Electron];
	dPIDTimingCuts[MuonPlus] = dPIDTimingCuts[Electron];

	// Timing Cuts: Mesons
	dPIDTimingCuts[PiPlus][SYS_BCAL] = 0.6;
	dPIDTimingCuts[PiPlus][SYS_FCAL] = 1.4;
	dPIDTimingCuts[PiPlus][SYS_TOF] = 0.3;

	dPIDTimingCuts[PiMinus][SYS_BCAL] = 0.6;
	dPIDTimingCuts[PiMinus][SYS_FCAL] = 1.4;
	dPIDTimingCuts[PiMinus][SYS_TOF] = 0.3;

	dPIDTimingCuts[KPlus][SYS_BCAL] = 0.4;
	dPIDTimingCuts[KPlus][SYS_FCAL] = 0.5;
	dPIDTimingCuts[KPlus][SYS_TOF] = 0.25;
	dPIDTimingCuts[KMinus] = dPIDTimingCuts[KPlus];

	// Timing Cuts: Baryons
	dPIDTimingCuts[Proton][SYS_BCAL] = 0.8;
	dPIDTimingCuts[Proton][SYS_FCAL] = 1.0;
	dPIDTimingCuts[Proton][SYS_TOF] = 0.5;

	dPIDTimingCuts[AntiProton] = dPIDTimingCuts[Proton];
}

map<Particle_t, pair<double, double> > DReaction_factory_trackeff_missing::Get_DecayingPIDs_InvariantMass(DReaction* locReaction)
{
	//pair: cut (hist) edge low/high
	map<Particle_t, pair<double, double> > locDecayingPIDs;
	for(size_t loc_i = 1; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		Particle_t locDecayingPID = locReactionStep->Get_InitialParticleID();
		bool locMissingMassFlag = locReaction->Check_IfMissingDecayProduct(loc_i);
		if(locMissingMassFlag)
			continue;

		auto locPIDIterator = dInvariantMassCuts.find(locDecayingPID);
		if(locPIDIterator == dInvariantMassCuts.end())
			continue;
		auto locCutPair = locPIDIterator->second;
		locDecayingPIDs[locDecayingPID] = locCutPair;
	}
	return locDecayingPIDs;
}

map<Particle_t, pair<int, deque<Particle_t> > > DReaction_factory_trackeff_missing::Get_DecayingPIDs_MissingMass(DReaction* locReaction)
{
	//int: decay step index
	map<Particle_t, pair<int, deque<Particle_t> > > locDecayingPIDs;
	for(size_t loc_i = 1; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		Particle_t locDecayingPID = locReactionStep->Get_InitialParticleID();
		bool locMissingMassFlag = locReaction->Check_IfMissingDecayProduct(loc_i);
		if(!locMissingMassFlag)
			continue;

		//missing mass cut
		auto locPIDIterator = dMissingMassCuts.find(locDecayingPID);
		if(locPIDIterator == dMissingMassCuts.end())
			continue;

		//get production step, other PIDs in that step
		deque<Particle_t> locMissingMassOffOfPIDs;
		int locDecayFromStep = locReaction->Get_InitialParticleDecayFromIndices(loc_i).first;
		locReaction->Get_ReactionStep(locDecayFromStep)->Get_FinalParticleIDs(locMissingMassOffOfPIDs);
		for(size_t loc_j = 0; loc_j < locMissingMassOffOfPIDs.size(); ++loc_j)
		{
			if(locMissingMassOffOfPIDs[loc_j] != locDecayingPID)
				continue;
			locMissingMassOffOfPIDs.erase(locMissingMassOffOfPIDs.begin() + loc_j);
			break;
		}

		locDecayingPIDs[locDecayingPID] = pair<int, deque<Particle_t> >(locDecayFromStep, locMissingMassOffOfPIDs);
	}
	return locDecayingPIDs;
}

void DReaction_factory_trackeff_missing::Add_PreComboCuts(DReaction* locReaction)
{
	// Highly Recommended: decaying particle mass cuts, applied during DParticleComboBlueprint construction
	map<Particle_t, pair<double, double> > locDecayingPIDs = Get_DecayingPIDs_InvariantMass(locReaction);
	for(auto& locMapPair : locDecayingPIDs)
		locReaction->Set_InvariantMassCut(locMapPair.first, locMapPair.second.first, locMapPair.second.second);

	map<Particle_t, pair<int, deque<Particle_t> > > locDecayingPIDs_Missing = Get_DecayingPIDs_MissingMass(locReaction);
	for(auto& locMapPair : locDecayingPIDs_Missing)
	{
		auto locCutPair = dMissingMassCuts[locMapPair.first];
		int locDecayFromStep = locMapPair.second.first;
		auto locMissingMassOffOfPIDs = locMapPair.second.second;

		//add the cut
		if(locCutPair.first >= 0.0)
			locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, locDecayFromStep, locMissingMassOffOfPIDs, false, locCutPair.first, locCutPair.second));
		else
			locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, locDecayFromStep, locMissingMassOffOfPIDs, false, locCutPair.first, locCutPair.second));
	}

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)

	//get cut pair
	pair<double, double> locCutPair;
	Particle_t locMissingPID;
	bool locMissingParticleFlag = locReaction->Get_MissingPID(locMissingPID); //false if none missing
	if(!locMissingParticleFlag)
		locCutPair = dMissingMassCuts[Unknown];
	else
		locCutPair = dMissingMassCuts[locMissingPID];

	if(locCutPair.first >= 0.0)
		locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMass(locReaction, false, locCutPair.first, locCutPair.second));
	else
		locReaction->Add_ComboPreSelectionAction(new DCutAction_MissingMassSquared(locReaction, false, locCutPair.first, locCutPair.second));
}

void DReaction_factory_trackeff_missing::Add_PIDActions(DReaction* locReaction)
{
	//Histogram before cuts
	locReaction->Add_ComboPreSelectionAction(new DHistogramAction_PID(locReaction));

	//Get, loop over detected PIDs in reaction
	deque<Particle_t> locDetectedPIDs;
	locReaction->Get_DetectedFinalPIDs(locDetectedPIDs);
	for(auto locPID : locDetectedPIDs)
	{
		if(dPIDTimingCuts.find(locPID) == dPIDTimingCuts.end())
			continue; //PID timing cut not defined!

		//Add timing cuts //false: measured data
		for(auto locSystemPair : dPIDTimingCuts[locPID])
			locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, locSystemPair.second, locPID, locSystemPair.first));

		//For kaon candidates, cut tracks that don't have a matching hit in a timing detector
			//Kaons are rare, cut ~necessary to reduce high backgrounds
		if((locPID == KPlus) || (locPID == KMinus))
			locReaction->Add_ComboPreSelectionAction(new DCutAction_NoPIDHit(locReaction, locPID));
	}

	//Loose dE/dx cuts
	locReaction->Add_ComboPreSelectionAction(new DCustomAction_dEdxCut_trackeff(locReaction, false)); //false: focus on keeping signal
}

void DReaction_factory_trackeff_missing::Add_MassHistograms(DReaction* locReaction, bool locUseKinFitResultsFlag, string locBaseUniqueName)
{
//	bool locConstrainMassFlag = locFSInfo->intermediateMassFits();
//	if(locUseKinFitResultsFlag && locConstrainMassFlag)
//		return;

	//missing mass
	if(!locUseKinFitResultsFlag)
	{
		//get cut pair
		pair<double, double> locCutPair;
		Particle_t locMissingPID;
		bool locMissingParticleFlag = locReaction->Get_MissingPID(locMissingPID); //false if none missing
		if(!locMissingParticleFlag)
			locCutPair = dMissingMassCuts[Unknown];
		else
			locCutPair = dMissingMassCuts[locMissingPID];

		//determine #bins
		int locNumBins = int((locCutPair.second - locCutPair.first)*1000.0 + 0.001);
		if(locNumBins < 200)
			locNumBins *= 5; //get close to 1000 bins
		if(locNumBins < 500)
			locNumBins *= 2; //get close to 1000 bins

		//add histogram action
		if(locCutPair.first >= 0.0)
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, locNumBins, locCutPair.first, locCutPair.second, locBaseUniqueName));
		else
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, locNumBins, locCutPair.first, locCutPair.second, locBaseUniqueName));
	}

	//invariant mass
	// Highly Recommended: Very loose decaying particle mass cuts, applied during DParticleComboBlueprint construction
	map<Particle_t, pair<double, double> > locDecayingPIDs = Get_DecayingPIDs_InvariantMass(locReaction);
	for(auto& locMapPair : locDecayingPIDs)
	{
		//determine #bins
		int locNumBins = int((locMapPair.second.second - locMapPair.second.first)*1000.0 + 0.001);
		if(locNumBins < 200)
			locNumBins *= 5; //get close to 1000 bins
		if(locNumBins < 500)
			locNumBins *= 2; //get close to 1000 bins

		//build name string
		Particle_t locPID = locMapPair.first;
		string locActionUniqueName = string(ParticleType(locPID)) + string("_") + locBaseUniqueName;

		//add histogram action
		locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, locPID, locUseKinFitResultsFlag,
			locNumBins, locMapPair.second.first, locMapPair.second.second, locActionUniqueName));
	}

	map<Particle_t, pair<int, deque<Particle_t> > > locDecayingPIDs_Missing = Get_DecayingPIDs_MissingMass(locReaction);
	for(auto& locMapPair : locDecayingPIDs_Missing)
	{
		auto locCutPair = dMissingMassCuts[locMapPair.first];
		int locDecayFromStep = locMapPair.second.first;
		auto locMissingMassOffOfPIDs = locMapPair.second.second;

		//determine #bins
		int locNumBins = int((locCutPair.second - locCutPair.first)*1000.0 + 0.001);
		if(locNumBins < 200)
			locNumBins *= 5; //get close to 1000 bins
		if(locNumBins < 500)
			locNumBins *= 2; //get close to 1000 bins

		//build name string
		Particle_t locPID = locMapPair.first;
		string locActionUniqueName = string(ParticleType(locPID)) + string("_") + locBaseUniqueName;

		//add the cut
		if(locCutPair.first >= 0.0)
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, locDecayFromStep, locMissingMassOffOfPIDs, locUseKinFitResultsFlag, locNumBins, locCutPair.first, locCutPair.second));
		else
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, locDecayFromStep, locMissingMassOffOfPIDs, locUseKinFitResultsFlag, locNumBins, locCutPair.first, locCutPair.second));
	}
}
