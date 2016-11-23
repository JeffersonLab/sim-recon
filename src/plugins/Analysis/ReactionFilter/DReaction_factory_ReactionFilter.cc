// $Id$
//
//    File: DReaction_factory_ReactionFilter.cc
// Created: Mon Nov 21 17:54:40 EST 2016
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//


#include "DReaction_factory_ReactionFilter.h"

//------------------
// init
//------------------
jerror_t DReaction_factory_ReactionFilter::init(void)
{
	Define_LooseCuts();

	//Get input reactions
	map<string, string> locParameterMap; //parameter key, value
	gPARMS->GetParameters(locParameterMap, "ReactionFilter"); //gets all parameters with this filter at the beginning of the key
	for(auto locParamPair : locParameterMap)
	{
		//hack so that don't get warning message about no default
		string locFullParamName = string("ReactionFilter") + locParamPair.first;
		string locFSValue;
		gPARMS->SetDefaultParameter(locFullParamName, locFSValue);
		dFSInfos.push_back(new FSInfo(locFSValue));
	}

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DReaction_factory_ReactionFilter::brun(JEventLoop* locEventLoop, int32_t locRunNumber)
{
	vector<double> locBeamPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
	dBeamBunchPeriod = locBeamPeriodVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DReaction_factory_ReactionFilter::evnt(JEventLoop* locEventLoop, uint64_t locEventNumber)
{
	// DOCUMENTATION:
	// ANALYSIS library: https://halldweb1.jlab.org/wiki/index.php/GlueX_Analysis_Software
	// DReaction factory: https://halldweb1.jlab.org/wiki/index.php/Analysis_DReaction

	// Loop over channels
	for(auto locFSInfo : dFSInfos)
	{
		/*************************************************** Reaction Definition **************************************************/

		//Create reaction //create with a unique name for each DReaction object. CANNOT (!) be "Thrown"
		DReaction* locReaction = new DReaction(locFSInfo->ReactionName());
		Create_FirstStep(locReaction, locFSInfo);
		Create_DecaySteps(locReaction, locFSInfo);

		/**************************************************** Control Settings ****************************************************/

		// Highly Recommended: Set EventStore skim query (use with "eventstore" source)
			// This will skip creating particle combos for events that aren't in the skims you list
			// Query should be comma-separated list of skims to boolean-AND together
		//locReaction->Set_EventStoreSkims("myskim1,myskim2,myskim3"); //boolean-AND of skims

		// Recommended: Type of kinematic fit to perform (default is d_NoFit)
			//fit types are of type DKinFitType, an enum defined in sim-recon/src/libraries/ANALYSIS/DReaction.h
			//Options: d_NoFit (default), d_P4Fit, d_VertexFit, d_P4AndVertexFit
			//P4 fits automatically constrain decaying particle masses, unless they are manually disabled
		locReaction->Set_KinFitType(d_P4AndVertexFit);

		// Highly Recommended: When generating particle combinations, reject all beam photons that match to a different RF bunch
		locReaction->Set_MaxPhotonRFDeltaT(1.5*dBeamBunchPeriod); // +/- 1 bunch for sideband subtraction

		// Highly Recommended: Cut on number of extra "good" tracks. "Good" tracks are ones that survive the "PreSelect" (or user custom) factory.
			// Important: Keep cut large: Can have many ghost and accidental tracks that look "good"
		if(!locFSInfo->inclusive())
		{
			int locMaxExtraGoodTracks = locFSInfo->missingN() ? 3 : 2;
			locReaction->Set_MaxExtraGoodTracks(locMaxExtraGoodTracks);
		}

		// Highly Recommended: Enable ROOT TTree output for this DReaction
		// string is file name (must end in ".root"!!): doen't need to be unique, feel free to change
		string locTreeFileName = string("tree_") + locFSInfo->ReactionName() + string(".root");
		bool locSaveUnusedFlag = locFSInfo->inclusive();
		locReaction->Enable_TTreeOutput(locTreeFileName, locSaveUnusedFlag); //true/false: do/don't save unused hypotheses

		/************************************************** Pre-Combo Custom Cuts *************************************************/

		Add_PreComboCuts(locReaction, locFSInfo);

		/**************************************************** Analysis Actions ****************************************************/

		// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
			//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
			//Pre-defined actions can be found in ANALYSIS/DHistogramActions_*.h and ANALYSIS/DCutActions.h
			//If a histogram action is repeated, it should be created with a unique name (string) to distinguish them

		// PID
		Add_PIDActions(locReaction);

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PreKinFit");

		// KINEMATIC FIT
		locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
		locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); //require kinematic fit converges

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PostKinFit");
		Add_MassHistograms(locReaction, locFSInfo, true, "PostKinFit_KinFit");

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, false)); //false: measured data

		_data.push_back(locReaction); //Register the DReaction with the factory
	}

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DReaction_factory_ReactionFilter::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReactionStepPool.size(); ++loc_i)
		delete dReactionStepPool[loc_i]; //cleanup memory
	return NOERROR;
}

/*********************************************************** CREATE REACTION STEPS ************************************************************/

void DReaction_factory_ReactionFilter::Create_FirstStep(DReaction* locReaction, FSInfo* locFSInfo)
{
	// the first reaction step
	DReactionStep* locReactionStep = new DReactionStep();
	locReactionStep->Set_InitialParticleID(Gamma);
	locReactionStep->Set_TargetParticleID(Proton);

	// primary particles
	vector<Particle_t> locPIDs = locFSInfo->PIDs();
	for (auto locPID : locPIDs)
		locReactionStep->Add_FinalParticleID(locPID);

	// add a missing nucleon
	if (locFSInfo->missingN() && locFSInfo->totalCharge() == 0)
		locReactionStep->Add_FinalParticleID(Proton, true);
	if (locFSInfo->missingN() && locFSInfo->totalCharge() == 1)
		locReactionStep->Add_FinalParticleID(Neutron, true);

	// add an unknown particle
	if (locFSInfo->inclusive())
		locReactionStep->Add_FinalParticleID(Unknown, true);

	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep);
}

void DReaction_factory_ReactionFilter::Create_DecaySteps(DReaction* locReaction, FSInfo* locFSInfo)
{
	// Keep track of DReactionSteps that have been created for decaying particles
	// These can be re-used between DReactions, allowing the analysis library to save memory for combo steps
		//bool: true for mass constrained by kinfit, false if not
	map<pair<Particle_t, bool>, set<DReactionStep*> > locDecayStepMap_Available = dDecayStepMap_All; 
	bool locConstrainMassFlag = locFSInfo->intermediateMassFits();

	vector<Particle_t> locPIDs = locFSInfo->PIDs();
	for(auto locPID : locPIDs)
	{
		//see if a decay step has been previously created for this PID
		pair<Particle_t, bool> locStepIDPair(locPID, locConstrainMassFlag);
		auto locStepMapIterator = locDecayStepMap_Available.find(locStepIDPair);
		if(locStepMapIterator == locDecayStepMap_Available.end())
		{
			Create_DecayStep(locReaction, locFSInfo, locPID); //nope, create it
			continue;
		}

		//see if there is an available decay step for this PID
		auto& locStepSet = locStepMapIterator->second;
		if(locStepSet.empty()) //nope, create it
			Create_DecayStep(locReaction, locFSInfo, locPID);
		else //re-use step
		{
			locReaction->Add_ReactionStep(*(locStepSet.begin()));
			locStepSet.erase(locStepSet.begin()); //no longer available
		}
	}
}

void DReaction_factory_ReactionFilter::Create_DecayStep(DReaction* locReaction, FSInfo* locFSInfo, Particle_t locPID)
{
	DReactionStep* locReactionStep = nullptr;
	if (locPID == Pi0)
	{
		locReactionStep = new DReactionStep();
		locReactionStep->Set_InitialParticleID(Pi0);
		locReactionStep->Add_FinalParticleID(Gamma);
		locReactionStep->Add_FinalParticleID(Gamma);
	}
	else if (locPID == Eta)
	{
		locReactionStep = new DReactionStep();
		locReactionStep->Set_InitialParticleID(Eta);
		locReactionStep->Add_FinalParticleID(Gamma);
		locReactionStep->Add_FinalParticleID(Gamma);
	}
	else if (locPID == Lambda)
	{
		locReactionStep = new DReactionStep();
		locReactionStep->Set_InitialParticleID(Lambda);
		locReactionStep->Add_FinalParticleID(Proton);
		locReactionStep->Add_FinalParticleID(PiMinus);
	}
	else if (locPID == AntiLambda)
	{
		locReactionStep = new DReactionStep();
		locReactionStep->Set_InitialParticleID(AntiLambda);
		locReactionStep->Add_FinalParticleID(AntiProton);
		locReactionStep->Add_FinalParticleID(PiPlus);
	}
	else if (locPID == KShort)
	{
		locReactionStep = new DReactionStep();
		locReactionStep->Set_InitialParticleID(KShort);
		locReactionStep->Add_FinalParticleID(PiPlus);
		locReactionStep->Add_FinalParticleID(PiMinus);
	}
	else
		return;

	bool locConstrainMassFlag = locFSInfo->intermediateMassFits();
	locReactionStep->Set_KinFitConstrainInitMassFlag(locConstrainMassFlag);

	locReaction->Add_ReactionStep(locReactionStep);
	dReactionStepPool.push_back(locReactionStep);

	pair<Particle_t, bool> locStepIDPair(locPID, locConstrainMassFlag);
	dDecayStepMap_All[locStepIDPair].insert(locReactionStep);
}

/************************************************************** ACTIONS AND CUTS **************************************************************/

void DReaction_factory_ReactionFilter::Define_LooseCuts(void)
{
	// Missing Mass Cuts
	dMissingMassCuts[Unknown] = pair<double, double>(-0.1, 0.1);
	dMissingMassCuts[Proton] = pair<double, double>(0.5, 1.4);
	dMissingMassCuts[Neutron] = pair<double, double>(0.5, 1.4);

	// Invariant Mass Cuts: Mesons
	dInvariantMassCuts[Pi0] = pair<double, double>(0.08, 0.19);
	dInvariantMassCuts[KShort] = pair<double, double>(0.3, 0.7);
	dInvariantMassCuts[Eta] = pair<double, double>(0.3, 0.8);
	dInvariantMassCuts[omega] = pair<double, double>(0.4, 1.2);
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
	dPIDTimingCuts[Gamma][SYS_BCAL] = 3.0;
	dPIDTimingCuts[Gamma][SYS_FCAL] = 2.5;

	// Timing Cuts: Leptons
	dPIDTimingCuts[Electron][SYS_BCAL] = 1.0;
	dPIDTimingCuts[Electron][SYS_FCAL] = 2.5;
	dPIDTimingCuts[Electron][SYS_TOF] = 2.5;

	dPIDTimingCuts[Positron][SYS_BCAL] = 1.0;
	dPIDTimingCuts[Positron][SYS_FCAL] = 2.5;
	dPIDTimingCuts[Positron][SYS_TOF] = 2.5;

	dPIDTimingCuts[MuonMinus][SYS_BCAL] = 1.0;
	dPIDTimingCuts[MuonMinus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[MuonMinus][SYS_TOF] = 2.5;

	dPIDTimingCuts[MuonPlus][SYS_BCAL] = 1.0;
	dPIDTimingCuts[MuonPlus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[MuonPlus][SYS_TOF] = 2.5;

	// Timing Cuts: Mesons
	dPIDTimingCuts[PiPlus][SYS_BCAL] = 2.0;
	dPIDTimingCuts[PiPlus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[PiPlus][SYS_TOF] = 2.5;

	dPIDTimingCuts[PiMinus][SYS_BCAL] = 2.0;
	dPIDTimingCuts[PiMinus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[PiMinus][SYS_TOF] = 2.5;

	dPIDTimingCuts[KPlus][SYS_BCAL] = 0.75;
	dPIDTimingCuts[KPlus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[KPlus][SYS_TOF] = 2.0;

	dPIDTimingCuts[KMinus][SYS_BCAL] = 0.75;
	dPIDTimingCuts[KMinus][SYS_FCAL] = 2.5;
	dPIDTimingCuts[KMinus][SYS_TOF] = 2.0;

	// Timing Cuts: Baryons
	dPIDTimingCuts[Proton][SYS_BCAL] = 2.5;
	dPIDTimingCuts[Proton][SYS_FCAL] = 2.5;
	dPIDTimingCuts[Proton][SYS_TOF] = 2.5;

	dPIDTimingCuts[AntiProton][SYS_BCAL] = 2.5;
	dPIDTimingCuts[AntiProton][SYS_FCAL] = 2.5;
	dPIDTimingCuts[AntiProton][SYS_TOF] = 2.5;
}

void DReaction_factory_ReactionFilter::Add_PreComboCuts(DReaction* locReaction, FSInfo* locFSInfo)
{
	// Highly Recommended: Very loose invariant mass cuts, applied during DParticleComboBlueprint construction
	vector<Particle_t> locPIDs = locFSInfo->PIDs();
	for(auto locPID : locPIDs)
	{
		auto locPIDIterator = dInvariantMassCuts.find(locPID);
		if(locPIDIterator == dInvariantMassCuts.end())
			continue;

		auto locCutPair = locPIDIterator->second;
		locReaction->Set_InvariantMassCut(locPID, locCutPair.first, locCutPair.second);
	}

	// Highly Recommended: Very loose DAnalysisAction cuts, applied just after creating the combination (before saving it)
	if(locFSInfo->exclusive())
	{
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
}

void DReaction_factory_ReactionFilter::Add_PIDActions(DReaction* locReaction)
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
	locReaction->Add_ComboPreSelectionAction(new DCustomAction_dEdxCut(locReaction, false)); //false: focus on keeping signal

	//Histogram after cuts
	locReaction->Add_ComboPreSelectionAction(new DHistogramAction_PID(locReaction, "PostPIDCuts"));
}

void DReaction_factory_ReactionFilter::Add_MassHistograms(DReaction* locReaction, FSInfo* locFSInfo, bool locUseKinFitResultsFlag, string locBaseUniqueName)
{
	bool locConstrainMassFlag = locFSInfo->intermediateMassFits();
	if(locUseKinFitResultsFlag && locConstrainMassFlag)
		return;

	//missing mass
	if(locFSInfo->exclusive() && !locUseKinFitResultsFlag)
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
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMass(locReaction, false, locNumBins, locCutPair.first, locCutPair.second));
		else
			locReaction->Add_AnalysisAction(new DHistogramAction_MissingMassSquared(locReaction, false, locNumBins, locCutPair.first, locCutPair.second));
	}

	//invariant mass
	vector<Particle_t> locPIDs = locFSInfo->PIDs();
	for(auto locPID : locPIDs)
	{
		auto locPIDIterator = dInvariantMassCuts.find(locPID);
		if(locPIDIterator == dInvariantMassCuts.end())
			continue;

		auto locCutPair = locPIDIterator->second;

		//determine #bins
		//range = 0.11
		//num bins = 110
		//range / #bins = good number
		//#bins = range*good_number
		int locNumBins = int((locCutPair.second - locCutPair.first)*1000.0 + 0.001);
		if(locNumBins < 200)
			locNumBins *= 5; //get close to 1000 bins
		if(locNumBins < 500)
			locNumBins *= 2; //get close to 1000 bins

		//build name string
		string locActionUniqueName = string(ParticleType(locPID)) + string("_") + locBaseUniqueName;

		//add histogram action
		locReaction->Add_AnalysisAction(new DHistogramAction_InvariantMass(locReaction, locPID, locUseKinFitResultsFlag, 
			locNumBins, locCutPair.first, locCutPair.second, locActionUniqueName));
	}
}

