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
	dSourceComboP4Handler = new DSourceComboP4Handler(nullptr, false);

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
		locReaction->Set_NumPlusMinusRFBunches(1); // +/- 1 bunch for sideband subtraction

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

		/**************************************************** Analysis Actions ****************************************************/

		// Recommended: Analysis actions automatically performed by the DAnalysisResults factories to histogram useful quantities.
			//These actions are executed sequentially, and are executed on each surviving (non-cut) particle combination 
			//Pre-defined actions can be found in ANALYSIS/DHistogramActions_*.h and ANALYSIS/DCutActions.h
			//If a histogram action is repeated, it should be created with a unique name (string) to distinguish them

		//PID
		locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PreKinFit");

		// KINEMATIC FIT
		locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
		//locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); //require kinematic fit converges

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PostKinFit");
		Add_MassHistograms(locReaction, locFSInfo, true, "PostKinFit_KinFit");

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //false: measured data

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
	auto locMissingPID = !locFSInfo->missingN() ? Unknown : ((locFSInfo->totalCharge() == 0) ? Proton : Neutron);
	auto locReactionStep = new DReactionStep(Gamma, Proton, locFSInfo->PIDs(), locMissingPID, locFSInfo->inclusive(), false);
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
	if(locPID == Pi0)
		locReactionStep = new DReactionStep(Pi0, {Gamma, Gamma});
	else if(locPID == Eta)
		locReactionStep = new DReactionStep(Eta, {Gamma, Gamma});
	else if(locPID == Lambda)
		locReactionStep = new DReactionStep(Lambda, {Proton, PiMinus});
	else if(locPID == AntiLambda)
		locReactionStep = new DReactionStep(AntiLambda, {AntiProton, PiPlus});
	else if(locPID == KShort)
		locReactionStep = new DReactionStep(KShort, {PiPlus, PiMinus});
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

void DReaction_factory_ReactionFilter::Add_MassHistograms(DReaction* locReaction, FSInfo* locFSInfo, bool locUseKinFitResultsFlag, string locBaseUniqueName)
{
	bool locConstrainMassFlag = locFSInfo->intermediateMassFits();
	if(locUseKinFitResultsFlag && locConstrainMassFlag)
		return;

	//invariant mass
	vector<Particle_t> locPIDs = locFSInfo->PIDs();
	set<Particle_t> locPIDsUsed;
	for(auto locPID : locPIDs)
	{
		pair<float, float> locCutPair;
		if(!dSourceComboP4Handler->Get_InvariantMassCuts(locPID, locCutPair))
			continue;
		if(locPIDsUsed.find(locPID) != locPIDsUsed.end())
			continue; //already done!

		//determine #bins
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

		locPIDsUsed.insert(locPID);
	}

	//missing mass
	if(!locFSInfo->exclusive() || locUseKinFitResultsFlag)
		return;

	//get cut pair
	auto locMissingPIDs = locReaction->Get_MissingPIDs();
	if(locMissingPIDs.size() >= 2)
		return;

	auto locMissingPID = locMissingPIDs.empty() ? Unknown : locMissingPIDs[0];
	pair<TF1*, TF1*> locFuncPair;
	if(!dSourceComboP4Handler->Get_MissingMassSquaredCuts(locMissingPID, locFuncPair))
		return;

	auto locCutPair = std::make_pair(locFuncPair.first->Eval(12.0), locFuncPair.second->Eval(12.0));

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

