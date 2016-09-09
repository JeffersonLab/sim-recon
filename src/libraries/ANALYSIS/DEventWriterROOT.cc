#include "DEventWriterROOT.h"

DEventWriterROOT::DEventWriterROOT(JEventLoop* locEventLoop)
{
	dInitNumThrownArraySize = 20;
	dInitNumBeamArraySize = 20;
	dInitNumTrackArraySize = 50;
	dInitNumNeutralArraySize = 15;
	dInitNumComboArraySize = 100;
	dThrownTreeInterface = NULL;

	//BEWARE: IF THIS IS CHANGED, CHANGE IN THE BLUEPRINT FACTORY AND THE ANALYSIS UTILITIES ALSO!!
	dTrackSelectionTag = "PreSelect";
	dShowerSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	locEventLoop->GetSingle(dAnalysisUtilities);

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	//CREATE & INITIALIZE ANALYSIS ACTIONS
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		if(!locReactions[loc_i]->Get_EnableTTreeOutputFlag())
			continue;

		dCutActionMap_ThrownTopology[locReactions[loc_i]] = new DCutAction_ThrownTopology(locReactions[loc_i], true);
		dCutActionMap_ThrownTopology[locReactions[loc_i]]->Initialize(locEventLoop);

		dCutActionMap_TrueCombo[locReactions[loc_i]] = new DCutAction_TrueCombo(locReactions[loc_i], 5.73303E-7, true); //+/- 5sigma
		dCutActionMap_TrueCombo[locReactions[loc_i]]->Initialize(locEventLoop);

		dCutActionMap_BDTSignalCombo[locReactions[loc_i]] = new DCutAction_BDTSignalCombo(locReactions[loc_i], 5.73303E-7, true, true); //+/- 5sigma
		dCutActionMap_BDTSignalCombo[locReactions[loc_i]]->Initialize(locEventLoop);
	}

	//CREATE TREES
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	//Get Target Center Z
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	dTargetCenterZ = 65.0;
	locGeometry->GetTargetZ(dTargetCenterZ);

	//CREATE TTREES
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		if(locReactions[loc_i]->Get_EnableTTreeOutputFlag())
			Create_DataTree(locReactions[loc_i], locEventLoop, !locMCThrowns.empty());
	}
}

DEventWriterROOT::~DEventWriterROOT(void)
{
	//Delete tree interface objects
	for(auto& locMapPair : dTreeInterfaceMap)
		delete locMapPair.second;
	if(dThrownTreeInterface != NULL)
		delete dThrownTreeInterface;

	for(auto& locMapPair : dTreeFillDataMap)
		delete locMapPair.second;

	//Delete actions
	for(auto& locMapPair : dCutActionMap_ThrownTopology)
		delete locMapPair.second;
	for(auto& locMapPair : dCutActionMap_TrueCombo)
		delete locMapPair.second;
	for(auto& locMapPair : dCutActionMap_BDTSignalCombo)
		delete locMapPair.second;
}

void DEventWriterROOT::Create_ThrownTree(JEventLoop* locEventLoop, string locOutputFileName) const
{
	dThrownTreeInterface = DTreeInterface::Create_DTreeInterface("Thrown_Tree", locOutputFileName);

	//TTREE BRANCHES
	DTreeBranchRegister locBranchRegister;

	//Get target PID
	const DMCReaction* locMCReaction = NULL;
	locEventLoop->GetSingle(locMCReaction);
	Particle_t locTargetPID = locMCReaction->target.PID();

	//setup target info
	Create_UserTargetInfo(locBranchRegister, locTargetPID);

	//create basic/misc. tree branches (run#, event#, etc.)
	locBranchRegister.Register_Single<UInt_t>("RunNumber");
	locBranchRegister.Register_Single<ULong64_t>("EventNumber");

	//Thrown Data
	Create_Branches_Thrown(locBranchRegister, true);

	//CUSTOM
	Create_CustomBranches_ThrownTree(locBranchRegister, locEventLoop);

	//CREATE BRANCHES
	dThrownTreeInterface->Create_Branches(locBranchRegister);
}

void DEventWriterROOT::Create_DataTree(const DReaction* locReaction, JEventLoop* locEventLoop, bool locIsMCDataFlag)
{
	string locReactionName = locReaction->Get_ReactionName();
	string locOutputFileName = locReaction->Get_TTreeOutputFileName();
	string locTreeName = locReactionName + string("_Tree");

	//create fill object
	dTreeFillDataMap[locReaction] = new DTreeFillData();

	//create tree interface
	DTreeInterface* locTreeInterface = DTreeInterface::Create_DTreeInterface(locTreeName, locOutputFileName);
	dTreeInterfaceMap[locReaction] = locTreeInterface;
	if(locTreeInterface->Get_BranchesCreatedFlag())
		return; //branches already created, then return

	//Branch register
	DTreeBranchRegister locBranchRegister;

	//fill maps
	TMap* locPositionToNameMap = Create_UserInfoMaps(locBranchRegister, locReaction);

/******************************************************************** Create Branches ********************************************************************/

	//create basic/misc. tree branches (run#, event#, etc.)
	locBranchRegister.Register_Single<UInt_t>("RunNumber");
	locBranchRegister.Register_Single<ULong64_t>("EventNumber");
	locBranchRegister.Register_Single<UInt_t>("L1TriggerBits");

	//create X4_Production
	locBranchRegister.Register_Single<TLorentzVector>("X4_Production");

	//create thrown branches
	if(locIsMCDataFlag)
	{
		Create_Branches_Thrown(locBranchRegister, false);
		locBranchRegister.Register_Single<Bool_t>("IsThrownTopology");
	}

	bool locBeamUsedFlag = (locReaction->Get_ReactionStep(0)->Get_TargetParticleID() != Unknown);

	//create branches for final-state particle hypotheses
	if(locBeamUsedFlag)
		Create_Branches_Beam(locBranchRegister, locIsMCDataFlag);
	Create_Branches_NeutralHypotheses(locBranchRegister, locIsMCDataFlag);
	Create_Branches_ChargedHypotheses(locBranchRegister, locIsMCDataFlag);

	//create branches for combos
	Create_Branches_Combo(locBranchRegister, locReaction, locIsMCDataFlag, locPositionToNameMap);

	//Custom branches
	Create_CustomBranches_DataTree(locBranchRegister, locEventLoop, locReaction, locIsMCDataFlag);

	//Create branches
	locTreeInterface->Create_Branches(locBranchRegister);
}

TMap* DEventWriterROOT::Create_UserInfoMaps(DTreeBranchRegister& locBranchRegister, const DReaction* locReaction) const
{
	//kinfit type
	DKinFitType locKinFitType = locReaction->Get_KinFitType();

	//create & add reaction identification maps
	TList* locUserInfo = locBranchRegister.Get_UserInfo();
	TMap* locNameToPIDMap = new TMap();
	locNameToPIDMap->SetName("NameToPIDMap");
	locUserInfo->Add(locNameToPIDMap);

	TMap* locNameToPositionMap = new TMap(); //not filled for target or initial particles
	locNameToPositionMap->SetName("NameToPositionMap");
	locUserInfo->Add(locNameToPositionMap);

	TMap* locPositionToNameMap = new TMap(); //not filled for target or initial particles
	locPositionToNameMap->SetName("PositionToNameMap");
	locUserInfo->Add(locPositionToNameMap);

	TMap* locPositionToPIDMap = new TMap();
	locPositionToPIDMap->SetName("PositionToPIDMap");
	locUserInfo->Add(locPositionToPIDMap);

	TMap* locDecayProductMap = new TMap(); //excludes resonances!!! //excludes intermediate decays (e.g. if xi- -> pi-, lambda -> pi-, pi-, p: will be xi- -> pi-, pi-, p and no lambda decay present)
	locDecayProductMap->SetName("DecayProductMap"); //parent name string -> tlist of decay product name strings
	locUserInfo->Add(locDecayProductMap);

	TMap* locMiscInfoMap = new TMap();
	locMiscInfoMap->SetName("MiscInfoMap");
	locUserInfo->Add(locMiscInfoMap);

	TList* locParticleNameList = new TList();
	locParticleNameList->SetName("ParticleNameList");
	locUserInfo->Add(locParticleNameList);

	//Set some misc info
	ostringstream locKinFitTypeStream;
	locKinFitTypeStream << locKinFitType;
	locMiscInfoMap->Add(new TObjString("KinFitType"), new TObjString(locKinFitTypeStream.str().c_str()));

	string HDDM_DATA_VERSION_STRING = "";
	if(gPARMS->Exists("REST:DATAVERSIONSTRING"))
		gPARMS->GetParameter("REST:DATAVERSIONSTRING", HDDM_DATA_VERSION_STRING);

	if(HDDM_DATA_VERSION_STRING != "")
		locMiscInfoMap->Add(new TObjString("REST:DATAVERSIONSTRING"), new TObjString(HDDM_DATA_VERSION_STRING.c_str()));

	string REST_JANA_CALIB_CONTEXT = "";
	if(gPARMS->Exists("REST:JANACALIBCONTEXT"))
		gPARMS->GetParameter("REST:JANACALIBCONTEXT", REST_JANA_CALIB_CONTEXT);
	if(REST_JANA_CALIB_CONTEXT == "")
		gPARMS->GetParameter("JANA_CALIB_CONTEXT", REST_JANA_CALIB_CONTEXT);
	if(REST_JANA_CALIB_CONTEXT != "")
		locMiscInfoMap->Add(new TObjString("REST:JANACALIBCONTEXT"), new TObjString(REST_JANA_CALIB_CONTEXT.c_str()));

	//find the # particles of each pid
	map<Particle_t, unsigned int> locParticleNumberMap;
	map<Particle_t, unsigned int> locDecayingParticleNumberMap;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j)) //missing particle
				continue;
			Particle_t locPID = locFinalParticleIDs[loc_j];

			//see if decays in a future step
			int locDecayStepIndex = locReaction->Get_DecayStepIndex(loc_i, loc_j);
			if(locDecayStepIndex >= 0) //decaying
			{
				if(locDecayingParticleNumberMap.find(locPID) == locDecayingParticleNumberMap.end())
					locDecayingParticleNumberMap[locPID] = 1;
				else
					++locDecayingParticleNumberMap[locPID];
			}
			else //detected, not decaying
			{
				if(locParticleNumberMap.find(locPID) == locParticleNumberMap.end())
					locParticleNumberMap[locPID] = 1;
				else
					++locParticleNumberMap[locPID];
			}
		}
	}

	//Create map objects
	map<Particle_t, unsigned int> locParticleNumberMap_Current, locDecayingParticleNumberMap_Current;
	Particle_t locTargetPID = Unknown;
	TObjString *locObjString_PID, *locObjString_Position, *locObjString_ParticleName;
	map<int, string> locDecayingParticleNames; //key is step index where they decay
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		{
			ostringstream locPIDStream, locPositionStream, locParticleNameStream;
			Particle_t locPID = locReactionStep->Get_InitialParticleID();
			locPIDStream << PDGtype(locPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());

			locPositionStream << loc_i << "_" << -1;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());

			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
			if((loc_i == 0) && ((locPID == Gamma) || (locPID == Electron) || (locPID == Positron)))
			{
				locParticleNameStream << "ComboBeam";
				locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
				locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
				locParticleNameList->AddLast(locObjString_ParticleName);
			}
			else //decaying particle
			{
				if(loc_i == 0)
					locParticleNameStream << "Decaying" << Convert_ToBranchName(ParticleType(locPID));
				else //name already created for final particle: use it
					locParticleNameStream << locDecayingParticleNames[loc_i];
				locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
				if(loc_i == 0) //in first step
				{
					locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
					locParticleNameList->AddLast(locObjString_ParticleName);
				}
			}
			locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);
			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
		}

		//target particle
		Particle_t locTempTargetPID = locReactionStep->Get_TargetParticleID();
		if((loc_i == 0) && (locTempTargetPID != Unknown))
		{
			locTargetPID = locTempTargetPID;

			ostringstream locPIDStream, locPositionStream, locParticleNameStream;
			locPIDStream << PDGtype(locTargetPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());

			locPositionStream << loc_i << "_" << -2;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());

			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);

			locParticleNameStream << "Target";
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());

			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);

			locParticleNameList->AddLast(locObjString_ParticleName);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			ostringstream locPIDStream, locPositionStream;
			Particle_t locPID = locFinalParticleIDs[loc_j];
			locPIDStream << PDGtype(locPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());

			locPositionStream << loc_i << "_" << loc_j;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());

			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j)) //missing particle
			{
				ostringstream locParticleNameStream;
				locParticleNameStream << "Missing" << Convert_ToBranchName(ParticleType(locPID));
				locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
				locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
				locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);
				locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
				string locPIDName = locParticleNameStream.str() + string("__PID");
				locMiscInfoMap->Add(new TObjString(locPIDName.c_str()), locObjString_PID);
				ostringstream locMassStream;
				locMassStream << ParticleMass(locPID);
				string locMassName = locParticleNameStream.str() + string("__Mass");
				locMiscInfoMap->Add(new TObjString(locMassName.c_str()), new TObjString(locMassStream.str().c_str()));
				locParticleNameList->AddLast(locObjString_ParticleName);
				continue;
			}

			//build name
			ostringstream locParticleNameStream;
			//see if decays in a future step
			int locDecayStepIndex = locReaction->Get_DecayStepIndex(loc_i, loc_j);
			if(locDecayStepIndex >= 0) //decays
			{
				if(locDecayingParticleNumberMap_Current.find(locPID) == locDecayingParticleNumberMap_Current.end())
					locDecayingParticleNumberMap_Current[locPID] = 1;
				else
					++locDecayingParticleNumberMap_Current[locPID];

				locParticleNameStream << "Decaying" << Convert_ToBranchName(ParticleType(locPID));
				if(locDecayingParticleNumberMap[locPID] > 1)
					locParticleNameStream << locDecayingParticleNumberMap_Current[locPID];
				locDecayingParticleNames[locDecayStepIndex] = locParticleNameStream.str();
			}
			else
			{
				if(locParticleNumberMap_Current.find(locPID) == locParticleNumberMap_Current.end())
					locParticleNumberMap_Current[locPID] = 1;
				else
					++locParticleNumberMap_Current[locPID];

				locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
				if(locParticleNumberMap[locPID] > 1)
					locParticleNameStream << locParticleNumberMap_Current[locPID];
			}

			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
			locParticleNameList->AddLast(locObjString_ParticleName);

			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
			locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			if(locDecayStepIndex >= 0)
			{
				ostringstream locMassStream;
				locMassStream << ParticleMass(locPID);
				string locMassName = locParticleNameStream.str() + string("__Mass");
				locMiscInfoMap->Add(new TObjString(locMassName.c_str()), new TObjString(locMassStream.str().c_str()));
			}
		}
	}

	//setup target info
	Create_UserTargetInfo(locBranchRegister, locTargetPID);

	//fill decay product map
	deque<size_t> locSavedSteps;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		Particle_t locPID = locReactionStep->Get_InitialParticleID();
		if(loc_i == 0)
			continue;
		if(!IsFixedMass(locPID))
			continue; //don't save resonance decays to the map

		//check to see if this decay has already been saved
		bool locStepAlreadySavedFlag = false;
		for(size_t loc_j = 0; loc_j < locSavedSteps.size(); ++loc_j)
		{
			if(locSavedSteps[loc_j] != loc_i)
				continue;
			locStepAlreadySavedFlag = true;
			break;
		}
		if(locStepAlreadySavedFlag)
			continue;

		//construct the name 
		ostringstream locParticleNameStream;
		locParticleNameStream << locDecayingParticleNames[loc_i];
		locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());

		TList* locDecayProductNames = NULL;
		Get_DecayProductNames(locReaction, loc_i, locPositionToNameMap, locDecayProductNames, locSavedSteps);
		locDecayProductMap->Add(locObjString_ParticleName, locDecayProductNames); //parent name string -> tobjarray of decay product name strings		
	}

	return locPositionToNameMap;
}

void DEventWriterROOT::Get_DecayProductNames(const DReaction* locReaction, size_t locReactionStepIndex, TMap* locPositionToNameMap, TList*& locDecayProductNames, deque<size_t>& locSavedSteps) const
{
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locReactionStepIndex);

	if(locDecayProductNames == NULL)
		locDecayProductNames = new TList();

	deque<Particle_t> locFinalParticleIDs;
	locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
	for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
	{
		//see if decays in a future step //save and continue if doesn't decay
		int locDecayStepIndex = locReaction->Get_DecayStepIndex(locReactionStepIndex, loc_j);
		if(locDecayStepIndex < 0)
		{
			ostringstream locPositionStream;
			locPositionStream << locReactionStepIndex << "_" << loc_j;
			locDecayProductNames->AddLast(locPositionToNameMap->GetValue(locPositionStream.str().c_str()));
			continue;
		}

		//add decay products
		Get_DecayProductNames(locReaction, locDecayStepIndex, locPositionToNameMap, locDecayProductNames, locSavedSteps);
	}

	locSavedSteps.push_back(locReactionStepIndex);
}

void DEventWriterROOT::Create_UserTargetInfo(DTreeBranchRegister& locBranchRegister, Particle_t locTargetPID) const
{
	TList* locUserInfo = locBranchRegister.Get_UserInfo();
	TMap* locMiscInfoMap = (TMap*)locUserInfo->FindObject("MiscInfoMap");
	if(locMiscInfoMap == NULL)
	{
		locMiscInfoMap = new TMap();
		locMiscInfoMap->SetName("MiscInfoMap");
		locUserInfo->Add(locMiscInfoMap);
	}

	//PID
	ostringstream locPIDStream;
	locPIDStream << PDGtype(locTargetPID);
	locMiscInfoMap->Add(new TObjString("Target__PID"), new TObjString(locPIDStream.str().c_str()));

	//Mass
	ostringstream locMassStream;
	locMassStream << ParticleMass(locTargetPID);
	locMiscInfoMap->Add(new TObjString("Target__Mass"), new TObjString(locMassStream.str().c_str()));

	//X, Y
	ostringstream locZeroStream;
	locZeroStream << 0.0;
	TObjString* locObjString_Zero = new TObjString(locZeroStream.str().c_str());
	locMiscInfoMap->Add(new TObjString("Target__CenterX"), locObjString_Zero);
	locMiscInfoMap->Add(new TObjString("Target__CenterY"), locObjString_Zero);

	//Z
	ostringstream locPositionStream;
	locPositionStream << dTargetCenterZ;
	TObjString* locObjString_Position = new TObjString(locPositionStream.str().c_str());
	locMiscInfoMap->Add(new TObjString("Target__CenterZ"), locObjString_Position);
}

void DEventWriterROOT::Create_Branches_Thrown(DTreeBranchRegister& locBranchRegister, bool locIsOnlyThrownFlag) const
{
	//BEAM
	locBranchRegister.Register_Single<Int_t>(Build_BranchName("ThrownBeam", "PID"));
	locBranchRegister.Register_Single<TLorentzVector>(Build_BranchName("ThrownBeam", "X4")); //reported at target center
	locBranchRegister.Register_Single<TLorentzVector>(Build_BranchName("ThrownBeam", "P4"));

	//EVENT-WIDE INFO
	locBranchRegister.Register_Single<ULong64_t>("NumPIDThrown_FinalState"); //19 digits
	locBranchRegister.Register_Single<ULong64_t>("PIDThrown_Decaying"); //19 digits
	locBranchRegister.Register_Single<Float_t>("MCWeight");

	//PRODUCTS
	Create_Branches_ThrownParticles(locBranchRegister, locIsOnlyThrownFlag);
}

void DEventWriterROOT::Create_Branches_ThrownParticles(DTreeBranchRegister& locBranchRegister, bool locIsOnlyThrownFlag) const
{
	string locParticleBranchName = "Thrown";

	string locArraySizeString = "NumThrown";
	locBranchRegister.Register_Single<UInt_t>(locArraySizeString);

	//IDENTIFIERS
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "ParentIndex"), locArraySizeString, dInitNumThrownArraySize);
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "PID"), locArraySizeString, dInitNumThrownArraySize);
	if(!locIsOnlyThrownFlag)
	{
		locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "MatchID"), locArraySizeString, dInitNumThrownArraySize);
		locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "MatchFOM"), locArraySizeString, dInitNumThrownArraySize);
	}

	//KINEMATICS: THROWN //at the production vertex
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4"), dInitNumThrownArraySize);
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4"), dInitNumThrownArraySize);
}

void DEventWriterROOT::Create_Branches_Beam(DTreeBranchRegister& locBranchRegister, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumBeam";
	locBranchRegister.Register_Single<UInt_t>(locArraySizeString);

	string locParticleBranchName = "Beam";

	//IDENTIFIER
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "PID"), locArraySizeString, dInitNumBeamArraySize);
	if(locIsMCDataFlag)
		locBranchRegister.Register_FundamentalArray<Bool_t>(Build_BranchName(locParticleBranchName, "IsGenerator"), locArraySizeString, dInitNumBeamArraySize);

	//KINEMATICS: MEASURED //at the production vertex
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), dInitNumThrownArraySize);
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), dInitNumThrownArraySize);
}

void DEventWriterROOT::Create_Branches_ChargedHypotheses(DTreeBranchRegister& locBranchRegister, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumChargedHypos";
	locBranchRegister.Register_Single<UInt_t>(locArraySizeString);

	string locParticleBranchName = "ChargedHypo";

	//IDENTIFIERS / MATCHING
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "TrackID"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "PID"), locArraySizeString, dInitNumTrackArraySize);
	if(locIsMCDataFlag)
		locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "ThrownIndex"), locArraySizeString, dInitNumTrackArraySize);

	//KINEMATICS //at the production vertex
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), dInitNumTrackArraySize);
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), dInitNumTrackArraySize);

	//TRACKING INFO
	locBranchRegister.Register_FundamentalArray<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Tracking"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Tracking"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_DCdEdx"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_DCdEdx"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_CDC"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_FDC"), locArraySizeString, dInitNumTrackArraySize);

	//TIMING INFO
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "HitTime"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "RFDeltaTVar"), locArraySizeString, dInitNumTrackArraySize);

	//PID QUALITY
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Timing"), locArraySizeString, dInitNumTrackArraySize);

	//HIT ENERGY
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_TOF"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_ST"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCAL"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCALPreshower"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_FCAL"), locArraySizeString, dInitNumTrackArraySize);

	//SHOWER MATCHING:
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaPhi"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaZ"), locArraySizeString, dInitNumTrackArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackFCAL_DOCA"), locArraySizeString, dInitNumTrackArraySize);
}

void DEventWriterROOT::Create_Branches_NeutralHypotheses(DTreeBranchRegister& locBranchRegister, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumNeutralHypos";
	string locParticleBranchName = "NeutralHypo";
	locBranchRegister.Register_Single<UInt_t>(locArraySizeString);

	//IDENTIFIERS / MATCHING
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "NeutralID"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "PID"), locArraySizeString, dInitNumNeutralArraySize);
	if(locIsMCDataFlag)
		locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "ThrownIndex"), locArraySizeString, dInitNumNeutralArraySize);

	//KINEMATICS //is combo-dependent: P4 defined by X4, X4 defined by other tracks
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), dInitNumNeutralArraySize);
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), dInitNumNeutralArraySize);

	//PID QUALITY
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Timing"), locArraySizeString, dInitNumNeutralArraySize);

	//SHOWER INFO
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Shower"), dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCAL"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCALPreshower"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Energy_FCAL"), locArraySizeString, dInitNumNeutralArraySize);

	//NEARBY TRACKS
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaPhi"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaZ"), locArraySizeString, dInitNumNeutralArraySize);
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "TrackFCAL_DOCA"), locArraySizeString, dInitNumNeutralArraySize);

	//PHOTON PID INFO
		//Computed using DVertex (best estimate of reaction vertex using all "good" tracks)
		//Can be used to compute timing chisq //is invalid for non-photons, so computed assuming photon PID
		//Variance of X4_Measured.T() - RFTime, regardless of which RF bunch is chosen. //RF bunch is combo-depende
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "PhotonRFDeltaTVar"), locArraySizeString, dInitNumNeutralArraySize);
}

void DEventWriterROOT::Create_Branches_Combo(DTreeBranchRegister& locBranchRegister, const DReaction* locReaction, bool locIsMCDataFlag, TMap* locPositionToNameMap) const
{
	string locNumComboString = "NumCombos";
	locBranchRegister.Register_Single<UInt_t>(locNumComboString);

	//kinfit type
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);

	//Is-cut
	locBranchRegister.Register_FundamentalArray<Bool_t>("IsComboCut", locNumComboString, dInitNumComboArraySize);

	//create combo-dependent, particle-independent branches
	if(locIsMCDataFlag)
	{
		locBranchRegister.Register_FundamentalArray<Bool_t>("IsTrueCombo", locNumComboString, dInitNumComboArraySize);
		locBranchRegister.Register_FundamentalArray<Bool_t>("IsBDTSignalCombo", locNumComboString, dInitNumComboArraySize);
	}

	locBranchRegister.Register_FundamentalArray<Float_t>("RFTime_Measured", locNumComboString, dInitNumComboArraySize);
	if(locKinFitFlag)
	{
		locBranchRegister.Register_FundamentalArray<Float_t>("ChiSq_KinFit", locNumComboString, dInitNumComboArraySize);
		locBranchRegister.Register_FundamentalArray<UInt_t>("NDF_KinFit", locNumComboString, dInitNumComboArraySize);
		if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
			locBranchRegister.Register_FundamentalArray<Float_t>("RFTime_KinFit", locNumComboString, dInitNumComboArraySize);
	}

	map<Particle_t, unsigned int> locParticleNumberMap_Current;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		Particle_t locInitialPID = locReactionStep->Get_InitialParticleID();
		//should check to make sure the beam particle isn't missing...
		if((loc_i == 0) && (locReactionStep->Get_InitialParticleID() != Unknown))
			Create_Branches_BeamComboParticle(locBranchRegister, locInitialPID, locKinFitType);
		else //decaying
		{
			//get the branch name
			ostringstream locPositionStream;
			locPositionStream << loc_i << "_-1";
			TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
			string locParticleBranchName = (const char*)(locObjString->GetString());

			if(IsFixedMass(locInitialPID) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
				locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), dInitNumComboArraySize);
			if((loc_i == 0) || IsDetachedVertex(locInitialPID))
				locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4"), dInitNumComboArraySize);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			int locDecayStepIndex = locReaction->Get_DecayStepIndex(loc_i, loc_j);
			if(locDecayStepIndex >= 0)
				continue; //decaying particle

			//get the branch name
			ostringstream locPositionStream;
			locPositionStream << loc_i << "_" << loc_j;
			TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
			string locParticleBranchName = (const char*)(locObjString->GetString());

			//missing particle
			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j))
			{
				// missing particle
				if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
					locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), dInitNumComboArraySize);
				continue;
			}

			Particle_t locPID = locFinalParticleIDs[loc_j];
			if(ParticleCharge(locPID) == 0)
				Create_Branches_ComboNeutral(locBranchRegister, locParticleBranchName, locKinFitType);
			else
				Create_Branches_ComboTrack(locBranchRegister, locParticleBranchName, locKinFitType);
		}
	}
}

void DEventWriterROOT::Create_Branches_BeamComboParticle(DTreeBranchRegister& locBranchRegister, Particle_t locBeamPID, DKinFitType locKinFitType) const
{
	string locParticleBranchName = "ComboBeam";
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "BeamIndex"), locArraySizeString, dInitNumComboArraySize);

	//KINEMATICS: KINFIT //at the interaction vertex
	if(locKinFitType != d_NoFit)
	{
		if(((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit)) || (ParticleCharge(locBeamPID) != 0))
			locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
			locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), dInitNumComboArraySize);
	}
}

void DEventWriterROOT::Create_Branches_ComboTrack(DTreeBranchRegister& locBranchRegister, string locParticleBranchName, DKinFitType locKinFitType) const
{
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "ChargedIndex"), locArraySizeString, dInitNumComboArraySize);

	//KINFIT INFO //at the production vertex
	if(locKinFitType != d_NoFit)
	{
		locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
		{
			locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), dInitNumComboArraySize);
			locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_KinFit"), locArraySizeString, dInitNumComboArraySize);
			locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_KinFit"), locArraySizeString, dInitNumComboArraySize);
		}
	}
}

void DEventWriterROOT::Create_Branches_ComboNeutral(DTreeBranchRegister& locBranchRegister, string locParticleBranchName, DKinFitType locKinFitType) const
{
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	locBranchRegister.Register_FundamentalArray<Int_t>(Build_BranchName(locParticleBranchName, "NeutralIndex"), locArraySizeString, dInitNumComboArraySize);

	//KINEMATICS //is combo-dependent: P4 defined by X4, X4 defined by other tracks
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), dInitNumComboArraySize);
	locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), dInitNumComboArraySize);

	//PID QUALITY
	locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_Measured"), locArraySizeString, dInitNumComboArraySize);
	if(locParticleBranchName.substr(0, 6) == "Photon")
		locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_Measured"), locArraySizeString, dInitNumComboArraySize);

	//KINFIT INFO //at the production vertex
	if(locKinFitType != d_NoFit)
	{
		locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
		{
			locBranchRegister.Register_ClonesArray<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), dInitNumComboArraySize);
			locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_KinFit"), locArraySizeString, dInitNumComboArraySize);
			if(locParticleBranchName.substr(0, 6) == "Photon")
				locBranchRegister.Register_FundamentalArray<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_KinFit"), locArraySizeString, dInitNumComboArraySize);
		}
	}
}

void DEventWriterROOT::Get_Reactions(jana::JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const
{
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	locReactions.clear();
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		if(string(locFactory->Tag()) == "Thrown")
			continue;

		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		locReactions.insert(locReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}
}

void DEventWriterROOT::Fill_ThrownTree(JEventLoop* locEventLoop) const
{
	vector<const DMCThrown*> locMCThrowns_FinalState;
	locEventLoop->Get(locMCThrowns_FinalState, "FinalState");

	vector<const DMCThrown*> locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");

	const DMCReaction* locMCReaction = NULL;
	locEventLoop->GetSingle(locMCReaction);

	ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
	Compute_ThrownPIDInfo(locMCThrowns_FinalState, locMCThrowns_Decaying, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

	vector<const DMCThrown*> locMCThrownsToSave;
	map<const DMCThrown*, unsigned int> locThrownIndexMap;
	Group_ThrownParticles(locMCThrowns_FinalState, locMCThrowns_Decaying, locMCThrownsToSave, locThrownIndexMap);

	//primary event info
	dThrownTreeFillData.Fill_Single<UInt_t>("RunNumber", locEventLoop->GetJEvent().GetRunNumber());
	dThrownTreeFillData.Fill_Single<ULong64_t>("EventNumber", locEventLoop->GetJEvent().GetEventNumber());

	//throwns
	Fill_ThrownInfo(&dThrownTreeFillData, locMCReaction, locMCThrownsToSave, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

	//Custom Branches
	Fill_CustomBranches_ThrownTree(&dThrownTreeFillData, locEventLoop, locMCReaction, locMCThrownsToSave);

	//FILL TTREE
	dThrownTreeInterface->Fill(dThrownTreeFillData);
}

void DEventWriterROOT::Fill_DataTrees(JEventLoop* locEventLoop, string locDReactionTag) const
{
	if(locDReactionTag == "Thrown")
	{
		cout << "WARNING: CANNOT FILL THROWN TREE WITH THIS FUNCTION." << endl;
		return;
	}

	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector);

	vector<const DReaction*> locReactionsWithTag;
	locEventLoop->Get(locReactionsWithTag, locDReactionTag.c_str());

	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		deque<const DParticleCombo*> locPassedParticleCombos;
		locAnalysisResultsVector[loc_i]->Get_PassedParticleCombos(locPassedParticleCombos);
		if(locPassedParticleCombos.empty())
			continue;
		const DReaction* locReaction = locAnalysisResultsVector[loc_i]->Get_Reaction();
		if(!locReaction->Get_EnableTTreeOutputFlag())
			continue;
		bool locReactionFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locReactionsWithTag.size(); ++loc_j)
		{
			if(locReactionsWithTag[loc_j] != locReaction)
				continue;
			locReactionFoundFlag = true;
			break;
		}
		if(!locReactionFoundFlag)
			continue; //reaction not from this factory, continue
		Fill_DataTree(locEventLoop, locReaction, locPassedParticleCombos);
	}
}

void DEventWriterROOT::Fill_DataTree(JEventLoop* locEventLoop, const DReaction* locReaction, deque<const DParticleCombo*>& locParticleCombos) const
{
	if(locReaction->Get_ReactionName() == "Thrown")
	{
		cout << "WARNING: CANNOT FILL THROWN TREE WITH THIS FUNCTION." << endl;
		return;
	}
	
	if(!locReaction->Get_EnableTTreeOutputFlag())
	{
		cout << "WARNING: ROOT TTREE OUTPUT NOT ENABLED FOR THIS DREACTION (" << locReaction->Get_ReactionName() << ")" << endl;
		return;
	}

	bool locSaveUnusedFlag = locReaction->Get_SaveUnusedFlag();

	//GET THROWN INFO
	vector<const DMCThrown*> locMCThrowns_FinalState;
	locEventLoop->Get(locMCThrowns_FinalState, "FinalState");

	vector<const DMCThrown*> locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

	vector<const DMCReaction*> locMCReactions;
	locEventLoop->Get(locMCReactions);
	const DMCReaction* locMCReaction = locMCReactions.empty() ? NULL : locMCReactions[0];

	//Pre-compute
	ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
	Compute_ThrownPIDInfo(locMCThrowns_FinalState, locMCThrowns_Decaying, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

	//Pre-compute
	vector<const DMCThrown*> locMCThrownsToSave;
	map<const DMCThrown*, unsigned int> locThrownIndexMap;
	Group_ThrownParticles(locMCThrowns_FinalState, locMCThrowns_Decaying, locMCThrownsToSave, locThrownIndexMap);

	//GET DETECTOR MATCHES
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//GET DVERTEX
	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	//GET TRIGGER
	const DTrigger* locTrigger = NULL;
	locEventLoop->GetSingle(locTrigger);

	//Check whether beam is used in the combo
	bool locBeamUsedFlag = (locReaction->Get_ReactionStep(0)->Get_TargetParticleID() != Unknown);

	//GET BEAM PHOTONS
		//however, only fill with beam particles that are in the combos
	set<const DBeamPhoton*> locBeamPhotonSet;
	vector<const DBeamPhoton*> locBeamPhotons;
	for(size_t loc_j = 0; loc_j < locParticleCombos.size(); ++loc_j)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombos[loc_j]->Get_ParticleComboStep(0);
		const DKinematicData* locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		if(locKinematicData == NULL)
			continue;
		const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locKinematicData);
		if(locBeamPhoton == NULL)
			continue;
		if(locBeamPhotonSet.find(locBeamPhoton) != locBeamPhotonSet.end())
			continue; //already included
		locBeamPhotonSet.insert(locBeamPhoton);
		locBeamPhotons.push_back(locBeamPhoton);
	}

	//create map of particles to array index: //beam now, will fill others in Get_*Hypotheses functions
		//used for pointing combo particles to the appropriate measured-particle array index
		//for hypos, they are the preselect versions if they exist, else the combo versions (e.g. PID not in REST)
	map<string, map<oid_t, int> > locObjectToArrayIndexMap;
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		locObjectToArrayIndexMap["DBeamPhoton"][locBeamPhotons[loc_i]->id] = loc_i;

	//Get PIDs need for reaction
	set<Particle_t> locReactionPIDs;
	//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles
	deque<Particle_t> locDetectedPIDs;
	locReaction->Get_DetectedFinalPIDs(locDetectedPIDs, 0, false); //q+
	for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
		locReactionPIDs.insert(locDetectedPIDs[loc_j]);

	//GET HYPOTHESES
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	if(locSaveUnusedFlag)
	{
		locChargedTrackHypotheses = Get_ChargedHypotheses(locEventLoop, locReactionPIDs, locObjectToArrayIndexMap["DChargedTrackHypothesis"]);
		locNeutralParticleHypotheses = Get_NeutralHypotheses(locEventLoop, locReactionPIDs, locObjectToArrayIndexMap["DNeutralParticleHypothesis"]);
	}
	else
	{
		locChargedTrackHypotheses = Get_ChargedHypotheses_Used(locEventLoop, locParticleCombos, locObjectToArrayIndexMap["DChargedTrackHypothesis"]);
		locNeutralParticleHypotheses = Get_NeutralHypotheses_Used(locEventLoop, locParticleCombos, locObjectToArrayIndexMap["DNeutralParticleHypothesis"]);
	}

	//EXECUTE ANALYSIS ACTIONS
	Bool_t locIsThrownTopologyFlag = kFALSE;
	vector<Bool_t> locIsTrueComboFlags;
	vector<Bool_t> locIsBDTSignalComboFlags;
	if(locMCReaction != NULL)
	{
		DCutAction_ThrownTopology* locThrownTopologyAction = dCutActionMap_ThrownTopology.find(locReaction)->second;
		locIsThrownTopologyFlag = (*locThrownTopologyAction)(locEventLoop, NULL); //combo not used/needed
		for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
		{
			DCutAction_TrueCombo* locTrueComboAction = dCutActionMap_TrueCombo.find(locReaction)->second;
			locIsTrueComboFlags.push_back((*locTrueComboAction)(locEventLoop, locParticleCombos[loc_i]));

			DCutAction_BDTSignalCombo* locBDTSignalComboAction = dCutActionMap_BDTSignalCombo.find(locReaction)->second;
			locIsBDTSignalComboFlags.push_back((*locBDTSignalComboAction)(locEventLoop, locParticleCombos[loc_i]));
		}
	}

	//Get tree fill data
	DTreeFillData* locTreeFillData = dTreeFillDataMap.find(locReaction)->second;

	/***************************************************** FILL TTREE DATA *****************************************************/

	//PRIMARY EVENT INFO
	locTreeFillData->Fill_Single<UInt_t>("RunNumber", locEventLoop->GetJEvent().GetRunNumber());
	locTreeFillData->Fill_Single<ULong64_t>("EventNumber", locEventLoop->GetJEvent().GetEventNumber());
	locTreeFillData->Fill_Single<UInt_t>("L1TriggerBits", locTrigger->Get_L1TriggerBits());

	//PRODUCTION X4
	DLorentzVector locProductionX4 = locVertex->dSpacetimeVertex;
	TLorentzVector locProductionTX4(locProductionX4.X(), locProductionX4.Y(), locProductionX4.Z(), locProductionX4.T());
	locTreeFillData->Fill_Single<TLorentzVector>("X4_Production", locProductionTX4);

	//THROWN INFORMATION
	if(locMCReaction != NULL)
	{
		Fill_ThrownInfo(locTreeFillData, locMCReaction, locMCThrownsToSave, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying, locMCThrownMatching, locObjectToArrayIndexMap);
		locTreeFillData->Fill_Single<Bool_t>("IsThrownTopology", locIsThrownTopologyFlag);
	}

	//INDEPENDENT BEAM PARTICLES
	if(locBeamUsedFlag)
	{
		//however, only fill with beam particles that are in the combos
		locTreeFillData->Fill_Single<UInt_t>("NumBeam", UInt_t(locBeamPhotons.size()));
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			Fill_BeamData(locTreeFillData, loc_i, locBeamPhotons[loc_i], locVertex, locMCThrownMatching);
	}

	//INDEPENDENT CHARGED TRACKS
	locTreeFillData->Fill_Single<UInt_t>("NumChargedHypos", UInt_t(locChargedTrackHypotheses.size()));
	for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); ++loc_i)
		Fill_ChargedHypo(locTreeFillData, loc_i, locChargedTrackHypotheses[loc_i], locMCThrownMatching, locThrownIndexMap, locDetectorMatches);

	//INDEPENDENT NEUTRAL PARTICLES
	locTreeFillData->Fill_Single<UInt_t>("NumNeutralHypos", UInt_t(locNeutralParticleHypotheses.size()));
	for(size_t loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); ++loc_i)
		Fill_NeutralHypo(locTreeFillData, loc_i, locNeutralParticleHypotheses[loc_i], locMCThrownMatching, locThrownIndexMap, locDetectorMatches);

	//COMBOS
	locTreeFillData->Fill_Single<UInt_t>("NumCombos", UInt_t(locParticleCombos.size()));
	for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
	{
		Fill_ComboData(locTreeFillData, locParticleCombos[loc_i], loc_i, locObjectToArrayIndexMap);
		if(locMCReaction != NULL)
		{
			locTreeFillData->Fill_Array<Bool_t>("IsTrueCombo", locIsTrueComboFlags[loc_i], loc_i);
			locTreeFillData->Fill_Array<Bool_t>("IsBDTSignalCombo", locIsTrueComboFlags[loc_i], loc_i);
		}
	}

	//CUSTOM
	Fill_CustomBranches_DataTree(locTreeFillData, locEventLoop, locMCReaction, locMCThrownsToSave, locMCThrownMatching, locDetectorMatches, locBeamPhotons, locChargedTrackHypotheses, locNeutralParticleHypotheses, locParticleCombos);

	//FILL
	DTreeInterface* locTreeInterface = dTreeInterfaceMap.find(locReaction)->second;
	locTreeInterface->Fill(*locTreeFillData);
}

vector<const DChargedTrackHypothesis*> DEventWriterROOT::Get_ChargedHypotheses(JEventLoop* locEventLoop, set<Particle_t> locReactionPIDs, map<oid_t, int>& locObjectToArrayIndexMap) const
{
	//Want to save all "PreSelect" hypotheses to the tree, plus hypos from PIDs that were reconstructed later
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	vector<const DChargedTrackHypothesis*> locComboChargedTrackHypotheses;
	locEventLoop->Get(locComboChargedTrackHypotheses, "Combo");

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Search for hypotheses with PIDs that were needed for the DReaction but that weren't reconstructed:
	map<const DChargedTrack*, map<Particle_t, const DChargedTrackHypothesis*> > locNewChargedHypothesesMap;

	//hypos not-findable in PreSelect, but not saved either //not saved: different RF time than main, or #-votes different
	set<const DChargedTrackHypothesis*> locUnsavedNonPreselectHypotheses;

	//loop
	for(size_t loc_i = 0; loc_i < locComboChargedTrackHypotheses.size(); ++loc_i)
	{
		Particle_t locPID = locComboChargedTrackHypotheses[loc_i]->PID();
		if(locReactionPIDs.find(locPID) == locReactionPIDs.end())
			continue; //PID not needed for this DReaction, don't bother

		//find Charged hypotheses not derived from the PreSelect factory
		const DChargedTrackHypothesis* locOrigChargedTrackHypothesis = NULL;
		locComboChargedTrackHypotheses[loc_i]->GetSingle(locOrigChargedTrackHypothesis);
		if(locOrigChargedTrackHypothesis != NULL)
			continue; //PID in REST file

		//OK, this hypo must at least be registered in locObjectToArrayIndexMap

		//Of these, choose the ones with the RF bunch time identical to the "main" one
		const DEventRFBunch* locComboEventRFBunch = NULL;
		locComboChargedTrackHypotheses[loc_i]->GetSingle(locComboEventRFBunch);
		if(fabs(locComboEventRFBunch->dTime - locEventRFBunch->dTime) > 0.01)
		{
			locUnsavedNonPreselectHypotheses.insert(locComboChargedTrackHypotheses[loc_i]);
			continue; //not the same one!
		}

		//Get original DChargedTrack
		const DChargedTrack* locOrigChargedTrack = NULL;
		locComboChargedTrackHypotheses[loc_i]->GetSingle(locOrigChargedTrack);

		//Of these hypos, choose only one: can be multiple (different #-votes): Make sure haven't saved one already
		map<Particle_t, const DChargedTrackHypothesis*>& locNewChargedHyposPIDMap = locNewChargedHypothesesMap[locOrigChargedTrack];
		if(locNewChargedHyposPIDMap.find(locPID) != locNewChargedHyposPIDMap.end())
		{
			locUnsavedNonPreselectHypotheses.insert(locComboChargedTrackHypotheses[loc_i]);
			continue; //already found for this PID //different #-votes
		}

		//unique! save!
		locNewChargedHyposPIDMap[locPID] = locComboChargedTrackHypotheses[loc_i];
	}

	//Build vector of combo-independent charged hypotheses to save
	vector<const DChargedTrackHypothesis*> locIndependentChargedTrackHypotheses;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const vector<const DChargedTrackHypothesis*>& locTrackHypoVector = locChargedTracks[loc_i]->dChargedTrackHypotheses;
		locIndependentChargedTrackHypotheses.insert(locIndependentChargedTrackHypotheses.end(), locTrackHypoVector.begin(), locTrackHypoVector.end());

		map<Particle_t, const DChargedTrackHypothesis*>& locNewChargedHyposPIDMap = locNewChargedHypothesesMap[locChargedTracks[loc_i]];
		map<Particle_t, const DChargedTrackHypothesis*>::iterator locMapIterator = locNewChargedHyposPIDMap.begin();
		for(; locMapIterator != locNewChargedHyposPIDMap.end(); ++locMapIterator)
			locIndependentChargedTrackHypotheses.push_back(locMapIterator->second);
	}

	//Register array indices for the main (directly-saved, measured) hypos
	for(size_t loc_i = 0; loc_i < locIndependentChargedTrackHypotheses.size(); ++loc_i)
		locObjectToArrayIndexMap[locIndependentChargedTrackHypotheses[loc_i]->id] = loc_i;

	//Register array indices for the not-directly-saved hypos (point to other saved ones)
	set<const DChargedTrackHypothesis*>::iterator locSetIterator = locUnsavedNonPreselectHypotheses.begin();
	for(; locSetIterator != locUnsavedNonPreselectHypotheses.end(); ++locSetIterator)
	{
		const DChargedTrack* locOrigChargedTrack = NULL;
		(*locSetIterator)->GetSingle(locOrigChargedTrack);
		oid_t locSourceID = locNewChargedHypothesesMap[locOrigChargedTrack][(*locSetIterator)->PID()]->id;
		locObjectToArrayIndexMap[(*locSetIterator)->id] = locObjectToArrayIndexMap[locSourceID];
	}

	return locIndependentChargedTrackHypotheses;
}

vector<const DChargedTrackHypothesis*> DEventWriterROOT::Get_ChargedHypotheses_Used(JEventLoop* locEventLoop, deque<const DParticleCombo*>& locParticleCombos, map<oid_t, int>& locObjectToArrayIndexMap) const
{
	set<const DChargedTrackHypothesis*> locFoundHypos;
	for(auto locCombo : locParticleCombos)
	{
		deque<const DKinematicData*> locChargedParticles;
		locCombo->Get_DetectedFinalChargedParticles_Measured(locChargedParticles);

		for(auto locParticle : locChargedParticles)
			locFoundHypos.insert(static_cast<const DChargedTrackHypothesis*>(locParticle));
	}

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Search for hypotheses with PIDs that were needed for the DReaction but that weren't reconstructed:
	map<const DChargedTrack*, map<Particle_t, const DChargedTrackHypothesis*> > locNewChargedHypothesesMap;

	//hypos not-findable in PreSelect, but not saved either //not saved: different RF time than main, or #-votes different
	set<const DChargedTrackHypothesis*> locUnsavedNonPreselectHypotheses;

	vector<const DChargedTrackHypothesis*> locIndependentChargedTrackHypotheses;
	set<const DChargedTrackHypothesis*> locSavedHypos;
	for(auto locChargedHypo : locFoundHypos)
	{
		Particle_t locPID = locChargedHypo->PID();
		const DChargedTrackHypothesis* locOrigChargedTrackHypothesis = NULL;
		locChargedHypo->GetSingle(locOrigChargedTrackHypothesis);

		//see if PID in REST file
		if(locOrigChargedTrackHypothesis != NULL)
		{
			//yes, save this (if not already)
			if(locSavedHypos.find(locOrigChargedTrackHypothesis) != locSavedHypos.end())
				continue;
			locIndependentChargedTrackHypotheses.push_back(locOrigChargedTrackHypothesis);
			locSavedHypos.insert(locOrigChargedTrackHypothesis);
			locObjectToArrayIndexMap[locOrigChargedTrackHypothesis->id] = locIndependentChargedTrackHypotheses.size() - 1;
			continue;
		}

		//Get original DChargedTrack
		const DChargedTrack* locOrigChargedTrack = NULL;
		locChargedHypo->GetSingle(locOrigChargedTrack);

		//Of these hypos, choose only one: can be multiple (different #-votes): Make sure haven't saved one already
		map<Particle_t, const DChargedTrackHypothesis*>& locNewChargedHyposPIDMap = locNewChargedHypothesesMap[locOrigChargedTrack];
		if(locNewChargedHyposPIDMap.find(locPID) != locNewChargedHyposPIDMap.end())
		{
			locUnsavedNonPreselectHypotheses.insert(locChargedHypo);
			continue; //already found for this PID //different #-votes
		}

		//unique! save!
		locNewChargedHyposPIDMap[locPID] = locChargedHypo;
		locIndependentChargedTrackHypotheses.push_back(locChargedHypo);
		locSavedHypos.insert(locChargedHypo);
		locObjectToArrayIndexMap[locChargedHypo->id] = locIndependentChargedTrackHypotheses.size() - 1;
	}

	//Register array indices for the not-directly-saved hypos (point to other saved ones)
	set<const DChargedTrackHypothesis*>::iterator locSetIterator = locUnsavedNonPreselectHypotheses.begin();
	for(; locSetIterator != locUnsavedNonPreselectHypotheses.end(); ++locSetIterator)
	{
		const DChargedTrack* locOrigChargedTrack = NULL;
		(*locSetIterator)->GetSingle(locOrigChargedTrack);
		oid_t locSourceID = locNewChargedHypothesesMap[locOrigChargedTrack][(*locSetIterator)->PID()]->id;
		locObjectToArrayIndexMap[(*locSetIterator)->id] = locObjectToArrayIndexMap[locSourceID];
	}

	return locIndependentChargedTrackHypotheses;
}

vector<const DNeutralParticleHypothesis*> DEventWriterROOT::Get_NeutralHypotheses(JEventLoop* locEventLoop, set<Particle_t> locReactionPIDs, map<oid_t, int>& locObjectToArrayIndexMap) const
{
	//Want to save all "PreSelect" hypotheses to the tree, plus hypos from PIDs that were reconstructed later
	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	vector<const DNeutralParticleHypothesis*> locComboNeutralParticleHypotheses;
	locEventLoop->Get(locComboNeutralParticleHypotheses, "Combo");

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Search for hypotheses with PIDs that were needed for the DReaction but that weren't reconstructed:
	map<const DNeutralParticle*, map<Particle_t, const DNeutralParticleHypothesis*> > locNewNeutralHypothesesMap;

	//hypos not-findable in PreSelect, but not saved either //not saved: different RF time than main, or #-votes different
	set<const DNeutralParticleHypothesis*> locUnsavedNonPreselectHypotheses;
	for(size_t loc_i = 0; loc_i < locComboNeutralParticleHypotheses.size(); ++loc_i)
	{
		Particle_t locPID = locComboNeutralParticleHypotheses[loc_i]->PID();
		if(locReactionPIDs.find(locPID) == locReactionPIDs.end())
			continue; //PID not needed for this DReaction, don't bother

		//find Neutral hypotheses not derived from the PreSelect factory
		const DNeutralParticleHypothesis* locOrigNeutralParticleHypothesis = NULL;
		locComboNeutralParticleHypotheses[loc_i]->GetSingle(locOrigNeutralParticleHypothesis);
		if(locOrigNeutralParticleHypothesis != NULL)
			continue; //PID in REST file

		//OK, this hypo must at least be registered in locObjectToArrayIndexMap

		//Of these, choose the ones with the RF bunch time identical to the "main" one
		const DEventRFBunch* locComboEventRFBunch = NULL;
		locComboNeutralParticleHypotheses[loc_i]->GetSingle(locComboEventRFBunch);
		if(fabs(locComboEventRFBunch->dTime - locEventRFBunch->dTime) > 0.01)
		{
			locUnsavedNonPreselectHypotheses.insert(locComboNeutralParticleHypotheses[loc_i]);
			continue; //not the same one!
		}

		//Get original DNeutralParticle
		const DNeutralParticle* locOrigNeutralParticle = NULL;
		locComboNeutralParticleHypotheses[loc_i]->GetSingle(locOrigNeutralParticle);

		//Of these hypos, choose only one: can be multiple (different #-votes): Make sure haven't saved one already
		map<Particle_t, const DNeutralParticleHypothesis*>& locNewNeutralHyposPIDMap = locNewNeutralHypothesesMap[locOrigNeutralParticle];
		if(locNewNeutralHyposPIDMap.find(locPID) != locNewNeutralHyposPIDMap.end())
		{
			locUnsavedNonPreselectHypotheses.insert(locComboNeutralParticleHypotheses[loc_i]);
			continue; //already found for this PID //different #-votes
		}

		//unique! save!
		locNewNeutralHyposPIDMap[locPID] = locComboNeutralParticleHypotheses[loc_i];
	}

	//Build vector of combo-independent Neutral hypotheses to save
	vector<const DNeutralParticleHypothesis*> locIndependentNeutralParticleHypotheses;
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		const vector<const DNeutralParticleHypothesis*>& locParticleHypoVector = locNeutralParticles[loc_i]->dNeutralParticleHypotheses;
		locIndependentNeutralParticleHypotheses.insert(locIndependentNeutralParticleHypotheses.end(), locParticleHypoVector.begin(), locParticleHypoVector.end());

		map<Particle_t, const DNeutralParticleHypothesis*>& locNewNeutralHyposPIDMap = locNewNeutralHypothesesMap[locNeutralParticles[loc_i]];
		map<Particle_t, const DNeutralParticleHypothesis*>::iterator locMapIterator = locNewNeutralHyposPIDMap.begin();
		for(; locMapIterator != locNewNeutralHyposPIDMap.end(); ++locMapIterator)
			locIndependentNeutralParticleHypotheses.push_back(locMapIterator->second);
	}

	//Register array indices for the main (directly-saved, measured) hypos
	for(size_t loc_i = 0; loc_i < locIndependentNeutralParticleHypotheses.size(); ++loc_i)
		locObjectToArrayIndexMap[locIndependentNeutralParticleHypotheses[loc_i]->id] = loc_i;

	//Register array indices for the not-directly-saved hypos (point to other saved ones)
	set<const DNeutralParticleHypothesis*>::iterator locSetIterator = locUnsavedNonPreselectHypotheses.begin();
	for(; locSetIterator != locUnsavedNonPreselectHypotheses.end(); ++locSetIterator)
	{
		const DNeutralParticle* locOrigNeutralParticle = NULL;
		(*locSetIterator)->GetSingle(locOrigNeutralParticle);
		oid_t locSourceID = locNewNeutralHypothesesMap[locOrigNeutralParticle][(*locSetIterator)->PID()]->id;
		locObjectToArrayIndexMap[(*locSetIterator)->id] = locObjectToArrayIndexMap[locSourceID];
	}

	return locIndependentNeutralParticleHypotheses;
}

vector<const DNeutralParticleHypothesis*> DEventWriterROOT::Get_NeutralHypotheses_Used(JEventLoop* locEventLoop, deque<const DParticleCombo*>& locParticleCombos, map<oid_t, int>& locObjectToArrayIndexMap) const
{
	set<const DNeutralParticleHypothesis*> locFoundHypos;
	for(auto locCombo : locParticleCombos)
	{
		deque<const DKinematicData*> locNeutralParticles;
		locCombo->Get_DetectedFinalNeutralParticles_Measured(locNeutralParticles);

		for(auto locParticle : locNeutralParticles)
			locFoundHypos.insert(static_cast<const DNeutralParticleHypothesis*>(locParticle));
	}

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	//Search for hypotheses with PIDs that were needed for the DReaction but that weren't reconstructed:
	map<const DNeutralParticle*, map<Particle_t, const DNeutralParticleHypothesis*> > locNewNeutralHypothesesMap;

	//hypos not-findable in PreSelect, but not saved either //not saved: different RF time than main, or #-votes different
	set<const DNeutralParticleHypothesis*> locUnsavedNonPreselectHypotheses;

	vector<const DNeutralParticleHypothesis*> locIndependentNeutralParticleHypotheses;
	set<const DNeutralParticleHypothesis*> locSavedHypos;
	for(auto locNeutralHypo : locFoundHypos)
	{
		Particle_t locPID = locNeutralHypo->PID();
		const DNeutralParticleHypothesis* locOrigNeutralParticleHypothesis = NULL;
		locNeutralHypo->GetSingle(locOrigNeutralParticleHypothesis);

		//see if PID in REST file
		if(locOrigNeutralParticleHypothesis != NULL)
		{
			//yes, save this (if not already)
			if(locSavedHypos.find(locOrigNeutralParticleHypothesis) != locSavedHypos.end())
				continue;
			locIndependentNeutralParticleHypotheses.push_back(locOrigNeutralParticleHypothesis);
			locSavedHypos.insert(locOrigNeutralParticleHypothesis);
			locObjectToArrayIndexMap[locOrigNeutralParticleHypothesis->id] = locIndependentNeutralParticleHypotheses.size() - 1;
			continue;
		}

		//Get original DNeutralParticle
		const DNeutralParticle* locOrigNeutralParticle = NULL;
		locNeutralHypo->GetSingle(locOrigNeutralParticle);

		//Of these hypos, choose only one: can be multiple (different #-votes): Make sure haven't saved one already
		map<Particle_t, const DNeutralParticleHypothesis*>& locNewNeutralHyposPIDMap = locNewNeutralHypothesesMap[locOrigNeutralParticle];
		if(locNewNeutralHyposPIDMap.find(locPID) != locNewNeutralHyposPIDMap.end())
		{
			locUnsavedNonPreselectHypotheses.insert(locNeutralHypo);
			continue; //already found for this PID //different #-votes
		}

		//unique! save!
		locNewNeutralHyposPIDMap[locPID] = locNeutralHypo;
		locIndependentNeutralParticleHypotheses.push_back(locNeutralHypo);
		locSavedHypos.insert(locNeutralHypo);
		locObjectToArrayIndexMap[locNeutralHypo->id] = locIndependentNeutralParticleHypotheses.size() - 1;
	}

	//Register array indices for the not-directly-saved hypos (point to other saved ones)
	set<const DNeutralParticleHypothesis*>::iterator locSetIterator = locUnsavedNonPreselectHypotheses.begin();
	for(; locSetIterator != locUnsavedNonPreselectHypotheses.end(); ++locSetIterator)
	{
		const DNeutralParticle* locOrigNeutralParticle = NULL;
		(*locSetIterator)->GetSingle(locOrigNeutralParticle);
		oid_t locSourceID = locNewNeutralHypothesesMap[locOrigNeutralParticle][(*locSetIterator)->PID()]->id;
		locObjectToArrayIndexMap[(*locSetIterator)->id] = locObjectToArrayIndexMap[locSourceID];
	}

	return locIndependentNeutralParticleHypotheses;
}

ULong64_t DEventWriterROOT::Calc_ParticleMultiplexID(Particle_t locPID) const
{
	int locPower = ParticleMultiplexPower(locPID);
	if(locPower == -1)
		return 0;

	int locIsFinalStateInt = Is_FinalStateParticle(locPID);
	if(locPID == Pi0)
		locIsFinalStateInt = 1;

	if(locIsFinalStateInt == 1) //decimal
	{
		ULong64_t locParticleMultiplexID = 1;
		for(int loc_i = 0; loc_i < locPower; ++loc_i)
			locParticleMultiplexID *= ULong64_t(10);
		return locParticleMultiplexID;
	}
	//decaying: binary
	return (ULong64_t(1) << ULong64_t(locPower)); //bit-shift
}

void DEventWriterROOT::Compute_ThrownPIDInfo(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying, ULong64_t& locNumPIDThrown_FinalState, ULong64_t& locPIDThrown_Decaying) const
{
	//THROWN PARTICLES BY PID
	locNumPIDThrown_FinalState = 0;
	for(size_t loc_i = 0; loc_i < locMCThrowns_FinalState.size(); ++loc_i) //final state
	{
		Particle_t locPID = locMCThrowns_FinalState[loc_i]->PID();
		ULong64_t locPIDMultiplexID = Calc_ParticleMultiplexID(locPID);
		if(locPIDMultiplexID == 0)
			continue; //unrecognized PID!!!
		unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % ULong64_t(10);
		if(locCurrentNumParticles != 9)
			locNumPIDThrown_FinalState += locPIDMultiplexID;
	}

	locPIDThrown_Decaying = 0;
	for(size_t loc_i = 0; loc_i < locMCThrowns_Decaying.size(); ++loc_i) //decaying
	{
		Particle_t locPID = locMCThrowns_Decaying[loc_i]->PID();
		ULong64_t locPIDMultiplexID = Calc_ParticleMultiplexID(locPID);
		if(locPIDMultiplexID == 0)
			continue; //unrecognized PID!!!
		if(locPID != Pi0)
			locPIDThrown_Decaying |= locPIDMultiplexID; //bit-wise or
		else //save pi0's as final state instead of decaying
		{
			unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % ULong64_t(10);
			if(locCurrentNumParticles != 9)
				locNumPIDThrown_FinalState += locPIDMultiplexID;
		}
	}
}

void DEventWriterROOT::Group_ThrownParticles(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying, vector<const DMCThrown*>& locMCThrownsToSave, map<const DMCThrown*, unsigned int>& locThrownIndexMap) const
{
	locMCThrownsToSave.clear();
	locMCThrownsToSave.insert(locMCThrownsToSave.end(), locMCThrowns_FinalState.begin(), locMCThrowns_FinalState.end());
	locMCThrownsToSave.insert(locMCThrownsToSave.end(), locMCThrowns_Decaying.begin(), locMCThrowns_Decaying.end());

	//create map of mcthrown to array index
	locThrownIndexMap.clear();
	for(size_t loc_i = 0; loc_i < locMCThrownsToSave.size(); ++loc_i)
		locThrownIndexMap[locMCThrownsToSave[loc_i]] = loc_i;
}

void DEventWriterROOT::Fill_ThrownInfo(DTreeFillData* locTreeFillData, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying) const
{
	map<string, map<oid_t, int> > locObjectToArrayIndexMap;
	Fill_ThrownInfo(locTreeFillData, locMCReaction, locMCThrowns, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying, NULL, locObjectToArrayIndexMap);
}

void DEventWriterROOT::Fill_ThrownInfo(DTreeFillData* locTreeFillData, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying, const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	//THIS MUST BE CALLED FROM WITHIN A LOCK, SO DO NOT PASS IN JEVENTLOOP! //TOO TEMPTING TO DO SOMETHING BAD

	//WEIGHT
	locTreeFillData->Fill_Single<Float_t>("MCWeight", locMCReaction->weight);

	//THROWN BEAM
	locTreeFillData->Fill_Single<Int_t>(Build_BranchName("ThrownBeam", "PID"), PDGtype(locMCReaction->beam.PID()));

	DVector3 locThrownBeamX3 = locMCReaction->beam.position();
	TLorentzVector locThrownBeamTX4(locThrownBeamX3.X(), locThrownBeamX3.Y(), locThrownBeamX3.Z(), locMCReaction->beam.time());
	locTreeFillData->Fill_Single<TLorentzVector>(Build_BranchName("ThrownBeam", "X4"), locThrownBeamTX4);

	DLorentzVector locThrownBeamP4 = locMCReaction->beam.lorentzMomentum();
	TLorentzVector locThrownBeamTP4(locThrownBeamP4.Px(), locThrownBeamP4.Py(), locThrownBeamP4.Pz(), locThrownBeamP4.E());
	locTreeFillData->Fill_Single<TLorentzVector>(Build_BranchName("ThrownBeam", "P4"), locThrownBeamTP4);

	//THROWN PRODUCTS
	locTreeFillData->Fill_Single<UInt_t>("NumThrown", locMCThrowns.size());
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		Fill_ThrownParticleData(locTreeFillData, loc_i, locMCThrowns[loc_i], locThrownIndexMap, locMCThrownMatching, locObjectToArrayIndexMap);

	//PID INFO
	locTreeFillData->Fill_Single<ULong64_t>("NumPIDThrown_FinalState", locNumPIDThrown_FinalState);
	locTreeFillData->Fill_Single<ULong64_t>("PIDThrown_Decaying", locPIDThrown_Decaying);
}

void DEventWriterROOT::Fill_ThrownParticleData(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DMCThrown* locMCThrown, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	string locParticleBranchName = "Thrown";

	//IDENTIFIERS
	int locParentIndex = -1; //e.g. photoproduced
	map<const DMCThrown*, unsigned int>::const_iterator locIterator;
	for(locIterator = locThrownIndexMap.begin(); locIterator != locThrownIndexMap.end(); ++locIterator)
	{
		if(locIterator->first->myid != locMCThrown->parentid)
			continue;
		locParentIndex = locIterator->second;
		break;
	}
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "ParentIndex"), locParentIndex, locArrayIndex);
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "PID"), locMCThrown->pdgtype, locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Int_t locMatchID = -1;
		double locMatchFOM = -1.0;
		if(ParticleCharge(locMCThrown->PID()) != 0)
		{
			const DChargedTrack* locChargedTrack = locMCThrownMatching->Get_MatchingChargedTrack(locMCThrown, locMatchFOM);
			if(locChargedTrack != NULL)
				locMatchID = locChargedTrack->Get_BestFOM()->candidateid;
		}
		else
		{
			//Can't use DNeutralShower JObject::id (must use dShowerID): 
				//Matching done with default-tag showers, but pre-select showers are saved to tree: JObject::id's don't match
			const DNeutralShower* locNeutralShower = locMCThrownMatching->Get_MatchingNeutralShower(locMCThrown, locMatchFOM);
			if(locNeutralShower != NULL)
				locMatchID = locNeutralShower->dShowerID;
		}
		locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "MatchID"), locMatchID, locArrayIndex);
		locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "MatchFOM"), locMatchFOM, locArrayIndex);
	}

	//KINEMATICS: THROWN //at the production vertex
	TLorentzVector locX4_Thrown(locMCThrown->position().X(), locMCThrown->position().Y(), locMCThrown->position().Z(), locMCThrown->time());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4"), locX4_Thrown, locArrayIndex);
	TLorentzVector locP4_Thrown(locMCThrown->momentum().X(), locMCThrown->momentum().Y(), locMCThrown->momentum().Z(), locMCThrown->energy());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4"), locP4_Thrown, locArrayIndex);
}

void DEventWriterROOT::Fill_BeamData(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DBeamPhoton* locBeamPhoton, const DVertex* locVertex, const DMCThrownMatching* locMCThrownMatching) const
{
	string locParticleBranchName = "Beam";

	//IDENTIFIER
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "PID"), PDGtype(locBeamPhoton->PID()), locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Bool_t locIsGeneratorFlag = (locMCThrownMatching->Get_ReconMCGENBeamPhoton() == locBeamPhoton) ? kTRUE : kFALSE;
		locTreeFillData->Fill_Array<Bool_t>(Build_BranchName(locParticleBranchName, "IsGenerator"), locIsGeneratorFlag, locArrayIndex);
	}

	//KINEMATICS: MEASURED

	//use production vertex, propagate photon time
	DVector3 locTargetCenter = locBeamPhoton->position();
	DVector3 locProductionVertex = locVertex->dSpacetimeVertex.Vect();
	double locDeltaPath = (locProductionVertex - locTargetCenter).Mag();
	bool locDownstreamFlag = ((locProductionVertex.Z() - locTargetCenter.Z()) > 0.0);
	double locDeltaT = locDownstreamFlag ? locDeltaPath/29.9792458 : -1.0*locDeltaPath/29.9792458;
	double locTime = locBeamPhoton->time() + locDeltaT;

	TLorentzVector locX4_Measured(locProductionVertex.X(), locProductionVertex.Y(), locProductionVertex.Z(), locTime);
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), locX4_Measured, locArrayIndex);

	DLorentzVector locDP4 = locBeamPhoton->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), locP4_Measured, locArrayIndex);
}

void DEventWriterROOT::Fill_ChargedHypo(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const
{
	string locParticleBranchName = "ChargedHypo";

	//ASSOCIATED OBJECTS
	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);

	const DBCALShower* locBCALShower = NULL;
	if(locChargedTrackHypothesis->Get_BCALShowerMatchParams() != NULL)
		locBCALShower = locChargedTrackHypothesis->Get_BCALShowerMatchParams()->dBCALShower;

	const DFCALShower* locFCALShower = NULL;
	if(locChargedTrackHypothesis->Get_FCALShowerMatchParams() != NULL)
		locFCALShower = locChargedTrackHypothesis->Get_FCALShowerMatchParams()->dFCALShower;

	//IDENTIFIERS
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "TrackID"), locChargedTrackHypothesis->candidateid, locArrayIndex);
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "PID"), PDGtype(locChargedTrackHypothesis->PID()), locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Int_t locThrownIndex = -1;
		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
		if(locMCThrown != NULL)
			locThrownIndex = locThrownIndexMap.find(locMCThrown)->second;
		locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "ThrownIndex"), locThrownIndex, locArrayIndex);
	}

	//KINEMATICS: MEASURED
	DVector3 locPosition = locChargedTrackHypothesis->position();
	TLorentzVector locTX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locChargedTrackHypothesis->time());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), locTX4_Measured, locArrayIndex);

	DLorentzVector locDP4 = locChargedTrackHypothesis->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), locP4_Measured, locArrayIndex);

	//TRACKING INFO
	locTreeFillData->Fill_Array<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Tracking"), locChargedTrackHypothesis->dNDF_Track, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Tracking"), locChargedTrackHypothesis->dChiSq_Track, locArrayIndex);
	locTreeFillData->Fill_Array<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_DCdEdx"), locChargedTrackHypothesis->dNDF_DCdEdx, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_DCdEdx"), locChargedTrackHypothesis->dChiSq_DCdEdx, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_CDC"), locTrackTimeBased->ddEdx_CDC, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_FDC"), locTrackTimeBased->ddEdx_FDC, locArrayIndex);

	//HIT ENERGY
	double locTOFdEdx = (locChargedTrackHypothesis->Get_TOFHitMatchParams() != NULL) ? locChargedTrackHypothesis->Get_TOFHitMatchParams()->dEdx : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_TOF"), locTOFdEdx, locArrayIndex);
	double locSCdEdx = (locChargedTrackHypothesis->Get_SCHitMatchParams() != NULL) ? locChargedTrackHypothesis->Get_SCHitMatchParams()->dEdx : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "dEdx_ST"), locSCdEdx, locArrayIndex);
	double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCAL"), locBCALEnergy, locArrayIndex);
	double locBCALPreshowerEnergy = (locBCALShower != NULL) ? locBCALShower->E_preshower : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCALPreshower"), locBCALPreshowerEnergy, locArrayIndex);

	double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_FCAL"), locFCALEnergy, locArrayIndex);

	//TIMING INFO
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "HitTime"), locChargedTrackHypothesis->t1(), locArrayIndex);
	double locStartTimeError = locChargedTrackHypothesis->t0_err();
	double locRFDeltaTVariance = (locChargedTrackHypothesis->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "RFDeltaTVar"), locRFDeltaTVariance, locArrayIndex);

	//MEASURED PID INFO
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing"), locChargedTrackHypothesis->measuredBeta(), locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing"), locChargedTrackHypothesis->dChiSq_Timing, locArrayIndex);
	locTreeFillData->Fill_Array<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Timing"), locChargedTrackHypothesis->dNDF_Timing, locArrayIndex);

	//SHOWER MATCHING: BCAL
	double locTrackBCAL_DeltaPhi = 999.0, locTrackBCAL_DeltaZ = 999.0;
	if(locChargedTrackHypothesis->Get_BCALShowerMatchParams() != NULL)
	{
		locTrackBCAL_DeltaPhi = locChargedTrackHypothesis->Get_BCALShowerMatchParams()->dDeltaPhiToShower;
		locTrackBCAL_DeltaZ = locChargedTrackHypothesis->Get_BCALShowerMatchParams()->dDeltaZToShower;
	}
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaPhi"), locTrackBCAL_DeltaPhi, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaZ"), locTrackBCAL_DeltaZ, locArrayIndex);

	//SHOWER MATCHING: FCAL
	double locDOCAToShower_FCAL = 999.0;
	if(locChargedTrackHypothesis->Get_FCALShowerMatchParams() != NULL)
		locDOCAToShower_FCAL = locChargedTrackHypothesis->Get_FCALShowerMatchParams()->dDOCAToShower;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackFCAL_DOCA"), locDOCAToShower_FCAL, locArrayIndex);
}

void DEventWriterROOT::Fill_NeutralHypo(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DMCThrownMatching* locMCThrownMatching, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const
{
	string locParticleBranchName = "NeutralHypo";
	const DNeutralShower* locNeutralShower = NULL;
	locNeutralParticleHypothesis->GetSingle(locNeutralShower);

	//ASSOCIATED OBJECTS
	const DBCALShower* locBCALShower = NULL;
	locNeutralShower->GetSingle(locBCALShower);
	const DFCALShower* locFCALShower = NULL;
	locNeutralShower->GetSingle(locFCALShower);

	//IDENTIFIERS
	Particle_t locPID = locNeutralParticleHypothesis->PID();
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "NeutralID"), locNeutralShower->dShowerID, locArrayIndex);
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "PID"), PDGtype(locPID), locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Int_t locThrownIndex = -1;
		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
		if(locMCThrown != NULL)
			locThrownIndex = locThrownIndexMap.find(locMCThrown)->second;
		locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "ThrownIndex"), locThrownIndex, locArrayIndex);
	}

	//KINEMATICS: MEASURED
	DVector3 locPosition = locNeutralParticleHypothesis->position();
	TLorentzVector locX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locNeutralParticleHypothesis->time());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), locX4_Measured, locArrayIndex);

	DLorentzVector locDP4 = locNeutralParticleHypothesis->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), locP4_Measured, locArrayIndex);

	//MEASURED PID INFO
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing"), locNeutralParticleHypothesis->measuredBeta(), locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing"), locNeutralParticleHypothesis->dChiSq, locArrayIndex);
	locTreeFillData->Fill_Array<UInt_t>(Build_BranchName(locParticleBranchName, "NDF_Timing"), locNeutralParticleHypothesis->dNDF, locArrayIndex);

	//SHOWER ENERGY
	DetectorSystem_t locDetector = locNeutralShower->dDetectorSystem;
	double locBCALEnergy = (locDetector == SYS_BCAL) ? locNeutralShower->dEnergy : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCAL"), locBCALEnergy, locArrayIndex);
	double locBCALPreshowerEnergy = (locDetector == SYS_BCAL) ? static_cast<const DBCALShower*>(locNeutralShower->dBCALFCALShower)->E_preshower : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_BCALPreshower"), locBCALPreshowerEnergy, locArrayIndex);
	double locFCALEnergy = (locDetector == SYS_FCAL) ? locNeutralShower->dEnergy : 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Energy_FCAL"), locFCALEnergy, locArrayIndex);

	//SHOWER POSITION
	DLorentzVector locHitDX4 = locNeutralShower->dSpacetimeVertex;
	TLorentzVector locTX4_Shower(locHitDX4.X(), locHitDX4.Y(), locHitDX4.Z(), locHitDX4.T());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Shower"), locTX4_Shower, locArrayIndex);

	//Track DOCA to Shower - BCAL
	double locNearestTrackBCALDeltaPhi = 999.0, locNearestTrackBCALDeltaZ = 999.0;
	if(locBCALShower != NULL)
	{
		if(!locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locNearestTrackBCALDeltaPhi, locNearestTrackBCALDeltaZ))
		{
			locNearestTrackBCALDeltaPhi = 999.0;
			locNearestTrackBCALDeltaZ = 999.0;
		}
		else if((locNearestTrackBCALDeltaPhi > 999.0) || (locNearestTrackBCALDeltaZ > 999.0))
		{
			locNearestTrackBCALDeltaPhi = 999.0;
			locNearestTrackBCALDeltaZ = 999.0;
		}
	}
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaPhi"), locNearestTrackBCALDeltaPhi, locArrayIndex);
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackBCAL_DeltaZ"), locNearestTrackBCALDeltaZ, locArrayIndex);

	//Track DOCA to Shower - FCAL
	double locDistanceToNearestTrack_FCAL = 999.0;
	if(locFCALShower != NULL)
	{
		if(!locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistanceToNearestTrack_FCAL))
			locDistanceToNearestTrack_FCAL = 999.0;
		if(locDistanceToNearestTrack_FCAL > 999.0)
			locDistanceToNearestTrack_FCAL = 999.0;
	}
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "TrackFCAL_DOCA"), locDistanceToNearestTrack_FCAL, locArrayIndex);

	//PHOTON PID INFO
	double locStartTimeError = locNeutralParticleHypothesis->t0_err();
	double locPhotonRFDeltaTVar = (locNeutralParticleHypothesis->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
	if(locPID != Gamma)
		locPhotonRFDeltaTVar = 0.0;
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "PhotonRFDeltaTVar"), locPhotonRFDeltaTVar, locArrayIndex);
}

void DEventWriterROOT::Fill_ComboData(DTreeFillData* locTreeFillData, const DParticleCombo* locParticleCombo, unsigned int locComboIndex, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	//MAIN CLASSES
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	//IS COMBO CUT
	locTreeFillData->Fill_Array<Bool_t>("IsComboCut", kFALSE, locComboIndex);

	//RF INFO
	double locRFTime = (locEventRFBunch != NULL) ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
	locTreeFillData->Fill_Array<Float_t>("RFTime_Measured", locRFTime, locComboIndex);

	//KINFIT INFO
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);
	if(locKinFitFlag)
	{
		if(locKinFitResults != NULL)
		{
			locTreeFillData->Fill_Array<Float_t>("ChiSq_KinFit", locKinFitResults->Get_ChiSq(), locComboIndex);
			locTreeFillData->Fill_Array<UInt_t>("NDF_KinFit", locKinFitResults->Get_NDF(), locComboIndex);
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
			{
				double locRFTime_KinFit = -9.9E9; //NOT IMPLEMENTED YET
				locTreeFillData->Fill_Array<Float_t>("RFTime_KinFit", locRFTime_KinFit, locComboIndex);
			}
		}
		else
		{
			locTreeFillData->Fill_Array<Float_t>("ChiSq_KinFit", 0.0, locComboIndex);
			locTreeFillData->Fill_Array<UInt_t>("NDF_KinFit", 0, locComboIndex);
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
				locTreeFillData->Fill_Array<Float_t>("RFTime_KinFit", -9.9E9, locComboIndex);
		}
	}

	//STEP DATA
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		Fill_ComboStepData(locTreeFillData, locParticleCombo, loc_i, locComboIndex, locKinFitType, locObjectToArrayIndexMap);
}

void DEventWriterROOT::Fill_ComboStepData(DTreeFillData* locTreeFillData, const DParticleCombo* locParticleCombo, unsigned int locStepIndex, unsigned int locComboIndex, DKinFitType locKinFitType, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	const TList* locUserInfo = dTreeInterfaceMap.find(locReaction)->second->Get_UserInfo();
	const TMap* locPositionToNameMap = (TMap*)locUserInfo->FindObject("PositionToNameMap");

	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	DLorentzVector locStepX4 = locParticleComboStep->Get_SpacetimeVertex();
	TLorentzVector locStepTX4(locStepX4.X(), locStepX4.Y(), locStepX4.Z(), locStepX4.T());

	//beam & production vertex
	Particle_t locInitialPID = locParticleComboStep->Get_InitialParticleID();
	const DKinematicData* locInitialParticle = locParticleComboStep->Get_InitialParticle();
	const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locInitialParticle);
	if(locBeamPhoton != NULL)
	{
		const DKinematicData* locInitParticleMeasured = locParticleComboStep->Get_InitialParticle_Measured();
		const DBeamPhoton* locMeasuredBeamPhoton = dynamic_cast<const DBeamPhoton*>(locInitParticleMeasured);
		int locBeamIndex = locObjectToArrayIndexMap.find("DBeamPhoton")->second.find(locMeasuredBeamPhoton->id)->second;
		Fill_ComboBeamData(locTreeFillData, locComboIndex, locBeamPhoton, locBeamIndex, locKinFitType);
	}
	else //decaying
	{
		//get the branch name
		ostringstream locPositionStream;
		locPositionStream << locStepIndex << "_-1";
		TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
		string locParticleBranchName = (const char*)(locObjString->GetString());

		if((locStepIndex == 0) || IsDetachedVertex(locInitialPID))
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4"), locStepTX4, locComboIndex);
		if(IsFixedMass(locInitialPID) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
		{
			TLorentzVector locDecayP4;
			if(locInitialParticle == NULL)
			{
				//fit failed to converge, calc from other particles
				DLorentzVector locDecayDP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, locStepIndex, false);
				locDecayDP4.SetE(sqrt(locDecayDP4.Vect().Mag2() + ParticleMass(locInitialPID)*ParticleMass(locInitialPID)));
				locDecayP4.SetPxPyPzE(locDecayDP4.Px(), locDecayDP4.Py(), locDecayDP4.Pz(), locDecayDP4.E());
			}
			else
				locDecayP4.SetPxPyPzE(locInitialParticle->momentum().X(), locInitialParticle->momentum().Y(), locInitialParticle->momentum().Z(), locInitialParticle->energy());
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), locDecayP4, locComboIndex);
		}
	}

	//final state particles
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_i);
		const DKinematicData* locKinematicData = locParticleComboStep->Get_FinalParticle(loc_i);
		const DKinematicData* locKinematicData_Measured = locParticleComboStep->Get_FinalParticle_Measured(loc_i);

		//decaying particle
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
			continue;

		//get the branch name
		ostringstream locPositionStream;
		locPositionStream << locStepIndex << "_" << loc_i;
		TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
		string locParticleBranchName = (const char*)(locObjString->GetString());

		//missing particle
		if(locParticleComboStep->Is_FinalParticleMissing(loc_i))
		{
			if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
			{
				TLorentzVector locMissingP4;
				if(locKinematicData == NULL)
				{
					//fit failed to converge, calc from other particles
					DLorentzVector locMissingDP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, false);
					locMissingDP4.SetE(sqrt(locMissingDP4.Vect().Mag2() + ParticleMass(locPID)*ParticleMass(locPID)));
					locMissingP4.SetPxPyPzE(locMissingDP4.Px(), locMissingDP4.Py(), locMissingDP4.Pz(), locMissingDP4.E());
				}
				else
					locMissingP4.SetPxPyPzE(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());

				locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), locMissingP4, locComboIndex);
			}
			continue;
		}

		//fill the data
		if(ParticleCharge(locPID) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData);
			const DNeutralParticleHypothesis* locMeasuredNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData_Measured);

			//find "NeutralIndex" //map contains PreSelect object if exists, else combo object
			const map<oid_t, int>& locObjectIDMap = locObjectToArrayIndexMap.find("DNeutralParticleHypothesis")->second;
			map<oid_t, int>::const_iterator locIDMapIterator = locObjectIDMap.find(locMeasuredNeutralHypo->id); //check for Combo
			if(locIDMapIterator == locObjectIDMap.end()) //not combo object: PreSelect
			{
				const DNeutralParticleHypothesis* locAssociatedNeutralHypo = NULL;
				locMeasuredNeutralHypo->GetSingle(locAssociatedNeutralHypo);
				locIDMapIterator = locObjectIDMap.find(locAssociatedNeutralHypo->id);
			}
			int locNeutralIndex = locIDMapIterator->second;

			Fill_ComboNeutralData(locTreeFillData, locComboIndex, locParticleBranchName, locMeasuredNeutralHypo, locNeutralHypo, locNeutralIndex, locKinFitType);
		}
		else
		{
			const DChargedTrackHypothesis* locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData);
			const DChargedTrackHypothesis* locMeasuredChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData_Measured);
			//find "ChargedIndex" //map contains PreSelect object if exists, else combo object
			const map<oid_t, int>& locObjectIDMap = locObjectToArrayIndexMap.find("DChargedTrackHypothesis")->second;
			map<oid_t, int>::const_iterator locIDMapIterator = locObjectIDMap.find(locMeasuredChargedHypo->id); //check for Combo
			if(locIDMapIterator == locObjectIDMap.end()) //not combo object: PreSelect
			{
				const DChargedTrackHypothesis* locAssociatedChargedHypo = NULL;
				locMeasuredChargedHypo->GetSingle(locAssociatedChargedHypo);
				locIDMapIterator = locObjectIDMap.find(locAssociatedChargedHypo->id);
			}
			int locChargedIndex = locIDMapIterator->second;

			Fill_ComboChargedData(locTreeFillData, locComboIndex, locParticleBranchName, locMeasuredChargedHypo, locChargedHypo, locChargedIndex, locKinFitType);
		}
	}
}

void DEventWriterROOT::Fill_ComboBeamData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, const DBeamPhoton* locBeamPhoton, unsigned int locBeamIndex, DKinFitType locKinFitType) const
{
	string locParticleBranchName = "ComboBeam";

	//IDENTIFIER
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "BeamIndex"), locBeamIndex, locComboIndex);

	//KINEMATICS: KINFIT
	if(locKinFitType != d_NoFit)
	{
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locBeamPhoton->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locBeamPhoton->time());
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), locX4_KinFit, locComboIndex);
		}

		//if charged, bends in b-field, update p4 when vertex changes
		if(((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit)) || (ParticleCharge(locBeamPhoton->PID()) != 0))
		{
			DLorentzVector locDP4 = locBeamPhoton->lorentzMomentum();
			TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), locP4_KinFit, locComboIndex);
		}
	}
}

void DEventWriterROOT::Fill_ComboChargedData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, string locParticleBranchName, const DChargedTrackHypothesis* locMeasuredChargedHypo, const DChargedTrackHypothesis* locChargedHypo, unsigned int locChargedIndex, DKinFitType locKinFitType) const
{
	//IDENTIFIER
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "ChargedIndex"), locChargedIndex, locComboIndex);

	//KINFIT
	if(locKinFitType != d_NoFit)
	{
		//KINEMATICS
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locChargedHypo->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locChargedHypo->time());
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), locX4_KinFit, locComboIndex);
		}

		//update even if vertex-only fit, because charged momentum propagated through b-field
		DLorentzVector locDP4 = locChargedHypo->lorentzMomentum();
		TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
		locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), locP4_KinFit, locComboIndex);

		//PID INFO
		if(locKinFitType != d_P4Fit)
		{
			locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_KinFit"), locChargedHypo->measuredBeta(), locComboIndex);
			locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_KinFit"), locChargedHypo->dChiSq_Timing, locComboIndex);
		}
	}
}

void DEventWriterROOT::Fill_ComboNeutralData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, string locParticleBranchName, const DNeutralParticleHypothesis* locMeasuredNeutralHypo, const DNeutralParticleHypothesis* locNeutralHypo, unsigned int locNeutralIndex, DKinFitType locKinFitType) const
{
	//IDENTIFIER
	locTreeFillData->Fill_Array<Int_t>(Build_BranchName(locParticleBranchName, "NeutralIndex"), locNeutralIndex, locComboIndex);

	//KINEMATICS: MEASURED
	DVector3 locPosition = locMeasuredNeutralHypo->position();
	TLorentzVector locX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locMeasuredNeutralHypo->time());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_Measured"), locX4_Measured, locComboIndex);

	DLorentzVector locDP4 = locMeasuredNeutralHypo->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_Measured"), locP4_Measured, locComboIndex);

	//MEASURED PID INFO
	locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_Measured"), locMeasuredNeutralHypo->measuredBeta(), locComboIndex);
	if(locParticleBranchName.substr(0, 6) == "Photon")
		locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_Measured"), locMeasuredNeutralHypo->dChiSq, locComboIndex);

	//KINFIT
	if(locKinFitType != d_NoFit)
	{
		//KINEMATICS
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locNeutralHypo->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locNeutralHypo->time());
			locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "X4_KinFit"), locX4_KinFit, locComboIndex);
		}

		//update even if vertex-only fit, because neutral momentum defined by vertex
		DLorentzVector locDP4 = locNeutralHypo->lorentzMomentum();
		TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
		locTreeFillData->Fill_Array<TLorentzVector>(Build_BranchName(locParticleBranchName, "P4_KinFit"), locP4_KinFit, locComboIndex);

		//PID INFO
		if(locKinFitType != d_P4Fit)
		{
			locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "Beta_Timing_KinFit"), locNeutralHypo->measuredBeta(), locComboIndex);
			if(locParticleBranchName.substr(0, 6) == "Photon")
				locTreeFillData->Fill_Array<Float_t>(Build_BranchName(locParticleBranchName, "ChiSq_Timing_KinFit"), locNeutralHypo->dChiSq, locComboIndex);
		}
	}
}
