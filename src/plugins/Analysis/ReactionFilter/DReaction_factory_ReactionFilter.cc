#include "DReaction_factory_ReactionFilter.h"

/*
 * COMMAND LINE INPUT:
 *
 * Any particles whose decays are not explicitly listed: default decay mode will be used if any. If none available, it will be assumed that the particle is to be detected.
 * This means that no default decay should be listed for the pi+, k+, etc.
 *
 * Reaction Form 1: BeamPID_TargetPID__FS1PID_FS2PID_FS3PID
 * Reaction Form 2: DecayPID_TargetPID__FS1PID_FS2PID_FS3PID
 * Reaction Form 3: DecayPID__FS1PID_FS2PID_FS3PID
 * To distinguish between forms 1 & 2: If first step, form 1 is assumed
 *
 * PID: Geant PID
 * parentheses around any means missing, Unknown means inclusive
 *
 * CONTROL MODE EXAMPLES:
 * -PReaction1=1_14__14_8_9 #g, p -> p, pi+, pi-
 * -PReaction2=1_14__14_8_9_7 #g, p -> p, pi+, pi-, pi0
 *
 * -PReaction3=1_14__14_8_9_7 #g, p -> p, pi+, pi-, pi0
 * -PReaction3:Decay1=7__2_3_1 # pi0 -> e+, e-, g
 * -PReaction3:Flags=F1_B2_T4_No7  #Fit enum 1 (p4-only (includes mass)), +/- 2 RF bunches, 4 extra tracks, don't constrain mass on pi0 (can have any number) # Can list these in any order
 *
 * -PReaction4=1_14__14_8_9_(7) #g, p -> p, pi+, pi-, missing pi0
 * -PReaction4=1_14__14_8_9_(0) #g, p -> p, pi+, pi-, inclusive
 *
 */

//------------------
// init
//------------------
jerror_t DReaction_factory_ReactionFilter::init(void)
{
	dSourceComboP4Handler = new DSourceComboP4Handler(nullptr, false);
	dSourceComboTimeHandler = new DSourceComboTimeHandler(nullptr, nullptr, nullptr);

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

/*
		//COMPARE:
        locReaction->Add_ComboPreSelectionAction(new DCutAction_NoPIDHit(locReaction, KPlus));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_NoPIDHit(locReaction, KMinus));

        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.5, Gamma, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Gamma, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Electron, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Electron, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Electron, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Positron, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Positron, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Positron, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiPlus, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiPlus, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiPlus, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.5, PiMinus, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, PiMinus, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, PiMinus, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, KPlus, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.75, KPlus, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, KPlus, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, KMinus, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.75, KMinus, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, KMinus, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.6, Proton, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, Proton, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, Proton, SYS_FCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 0.6, AntiProton, SYS_TOF));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 1.0, AntiProton, SYS_BCAL));
        locReaction->Add_ComboPreSelectionAction(new DCutAction_PIDDeltaT(locReaction, false, 2.0, AntiProton, SYS_FCAL));
*/

		//PID
		locReaction->Add_AnalysisAction(new DHistogramAction_PID(locReaction));

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PreKinFit");

		// KINEMATIC FIT
		locReaction->Add_AnalysisAction(new DHistogramAction_KinFitResults(locReaction, 0.05)); //5% confidence level cut on pull histograms only
		//locReaction->Add_AnalysisAction(new DCutAction_KinFitFOM(locReaction, -1.0)); //require kinematic fit converges

//bool Get_TimeCut(Particle_t locPID, DetectorSystem_t locSystem, TF1* locTimeCut_ns) const;

		// HISTOGRAM MASSES //false/true: measured/kinfit data
		Add_MassHistograms(locReaction, locFSInfo, false, "PostKinFit");
		Add_MassHistograms(locReaction, locFSInfo, true, "PostKinFit_KinFit");

		// KINEMATICS
		locReaction->Add_AnalysisAction(new DHistogramAction_ParticleComboKinematics(locReaction, true)); //true: kinfit data

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

	//MESONS
	if(locPID == Pi0)
		locReactionStep = new DReactionStep(Pi0, {Gamma, Gamma});
	else if(locPID == KShort)
		locReactionStep = new DReactionStep(KShort, {PiMinus, PiPlus});
	else if(locPID == Eta)
		locReactionStep = new DReactionStep(Eta, {Gamma, Gamma});
	else if(locPID == omega)
		locReactionStep = new DReactionStep(omega, {PiMinus, PiPlus, Pi0});
	else if(locPID == EtaPrime)
		locReactionStep = new DReactionStep(EtaPrime, {PiMinus, PiPlus, Eta});
	else if(locPID == phiMeson)
		locReactionStep = new DReactionStep(phiMeson, {KMinus, KPlus});
	else if(locPID == D0)
		locReactionStep = new DReactionStep(D0, {PiMinus, KPlus});
	else if(locPID == Jpsi)
		locReactionStep = new DReactionStep(Jpsi, {Electron, Positron});

	//BARYONS
	else if(locPID == Lambda)
		locReactionStep = new DReactionStep(Lambda, {PiMinus, Proton});
	else if(locPID == AntiLambda)
		locReactionStep = new DReactionStep(AntiLambda, {PiPlus, AntiProton});
	else if(locPID == SigmaMinus)
		locReactionStep = new DReactionStep(SigmaMinus, {PiMinus, Neutron});
	else if(locPID == Sigma0)
		locReactionStep = new DReactionStep(Sigma0, {Gamma, Lambda});
	else if(locPID == SigmaPlus)
		locReactionStep = new DReactionStep(SigmaPlus, {Pi0, Proton});
	else if(locPID == Xi0)
		locReactionStep = new DReactionStep(Xi0, {Pi0, Lambda});
	else if(locPID == XiMinus)
		locReactionStep = new DReactionStep(XiMinus, {PiMinus, Lambda});
	else if(locPID == OmegaMinus)
		locReactionStep = new DReactionStep(OmegaMinus, {KMinus, Lambda});
	else if(locPID == Lambda_c)
		locReactionStep = new DReactionStep(Lambda_c, {PiPlus, KMinus, Proton});
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

