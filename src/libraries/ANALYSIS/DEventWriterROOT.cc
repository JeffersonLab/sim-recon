#include "DEventWriterROOT.h"

DEventWriterROOT::DEventWriterROOT(JEventLoop* locEventLoop)
{
	dInitNumThrownArraySize = 20;
	dInitNumBeamArraySize = 20;
	dInitNumTrackArraySize = 50;
	dInitNumShowerArraySize = 15;
	dInitNumComboArraySize = 100;

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

	japp->RootWriteLock();
	{
		++Get_NumEventWriterThreads();
	}
	japp->RootUnLock();
}

DEventWriterROOT::~DEventWriterROOT(void)
{
	japp->RootWriteLock();
	{
		--Get_NumEventWriterThreads();
		if(Get_NumEventWriterThreads() != 0)
		{
			japp->RootUnLock();
			return;
		}
		for(size_t loc_i = 0; loc_i < Get_OutputROOTFiles().size(); ++loc_i)
		{
			Get_OutputROOTFiles()[loc_i]->Write(0, TObject::kOverwrite);
			Get_OutputROOTFiles()[loc_i]->Close();
			delete Get_OutputROOTFiles()[loc_i];
		}
		Get_OutputROOTFiles().clear();
	}
	japp->RootUnLock();
}

int& DEventWriterROOT::Get_NumEventWriterThreads(void) const
{
	static int locNumEventWriterThreads = 0;
	return locNumEventWriterThreads;
}

string& DEventWriterROOT::Get_ThrownTreeFileName(void) const
{
	//static so that it's not a member: can be changed in a call to a const function //object is const when the user gets it
	static string locThrownTreeFileName = "";
	return locThrownTreeFileName;
}

map<TTree*, map<string, TClonesArray*> >& DEventWriterROOT::Get_ClonesArrayMap(void) const
{
	static map<TTree*, map<string, TClonesArray*> > locClonesArrayMap;
	return locClonesArrayMap;
}

map<TTree*, map<string, TObject*> >& DEventWriterROOT::Get_TObjectMap(void) const
{
	static map<TTree*, map<string, TObject*> > locTObjectMap;
	return locTObjectMap;
}

map<TTree*, map<string, unsigned int> >& DEventWriterROOT::Get_FundamentalArraySizeMap(void) const
{
	static map<TTree*, map<string, unsigned int> > locFundamentalArraySizeMap;
	return locFundamentalArraySizeMap;
}

deque<TFile*>& DEventWriterROOT::Get_OutputROOTFiles(void) const
{
	static deque<TFile*> locOutputROOTFiles;
	return locOutputROOTFiles;
}

void DEventWriterROOT::Create_ThrownTree(string locOutputFileName) const
{
	japp->RootWriteLock();
	{
		if(Get_ThrownTreeFileName() != "")
		{
			japp->RootUnLock();
			return; //already created
		}
		Get_ThrownTreeFileName() = locOutputFileName;

		//create ROOT file if it doesn't exist already
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile == NULL)
		{
			locFile = new TFile(locOutputFileName.c_str(), "RECREATE");
			Get_OutputROOTFiles().push_back(locFile);
		}

		//create tree if it doesn't exist already
		string locTreeName = "Thrown_Tree";
		locFile->cd();
		if(gDirectory->Get(locTreeName.c_str()) != NULL)
		{
			japp->RootUnLock();
			return; //already created by another thread
		}
		TTree* locTree = new TTree(locTreeName.c_str(), locTreeName.c_str());

		/******************************************************************** Create Branches ********************************************************************/

		//create basic/misc. tree branches (run#, event#, etc.)
		Create_Branch_Fundamental<UInt_t>(locTree, "RunNumber");
		Create_Branch_Fundamental<ULong64_t>(locTree, "EventNumber");

		//Thrown Data
		Create_Branches_Thrown(locTree, true);

		//CUSTOM
		Create_CustomBranches_ThrownTree(locTree);
	}
	japp->RootUnLock();
}

void DEventWriterROOT::Create_DataTrees(JEventLoop* locEventLoop) const
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	//CREATE TTREES
	japp->RootWriteLock();
	{
		for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		{
			if(locReactions[loc_i]->Get_EnableTTreeOutputFlag())
				Create_DataTree(locReactions[loc_i], !locMCThrowns.empty());
		}
	}
	japp->RootUnLock();
}

void DEventWriterROOT::Create_DataTree(const DReaction* locReaction, bool locIsMCDataFlag) const
{
	string locReactionName = locReaction->Get_ReactionName();

	//create ROOT file if it doesn't exist already
	string locOutputFileName = locReaction->Get_TTreeOutputFileName();
	TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
	if(locFile == NULL)
	{
		locFile = new TFile(locOutputFileName.c_str(), "RECREATE");
		Get_OutputROOTFiles().push_back(locFile);
	}

	//create tree if it doesn't exist already
	string locTreeName = locReactionName + string("_Tree");
	locFile->cd();
	if(gDirectory->Get(locTreeName.c_str()) != NULL)
		return; //already created by another thread
	TTree* locTree = new TTree(locTreeName.c_str(), locTreeName.c_str());

	//find the # particles of each pid
	map<Particle_t, unsigned int> locParticleNumberMap;
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
			if(locParticleNumberMap.find(locPID) == locParticleNumberMap.end())
				locParticleNumberMap[locPID] = 1;
			else
				++locParticleNumberMap[locPID];
		}
	}

	//fill maps
	Create_UserInfoMaps(locTree, locReaction, locParticleNumberMap);

/******************************************************************** Create Branches ********************************************************************/

	//create basic/misc. tree branches (run#, event#, etc.)
	Create_Branch_Fundamental<UInt_t>(locTree, "RunNumber");
	Create_Branch_Fundamental<ULong64_t>(locTree, "EventNumber");

	//create thrown branches
	if(locIsMCDataFlag)
	{
		Create_Branches_Thrown(locTree, false);
		Create_Branch_Fundamental<Float_t>(locTree, "IsThrownTopology");
	}

	Particle_t locInitialPID = locReaction->Get_ReactionStep(0)->Get_InitialParticleID();
	bool locBeamUsedFlag = ((locInitialPID == Gamma) || (locInitialPID == Electron) || (locInitialPID == Positron));

	//create branches for final-state particle hypotheses
	if(locBeamUsedFlag)
		Create_Branches_Beam(locTree, locIsMCDataFlag);
	Create_Branches_NeutralShowers(locTree, locIsMCDataFlag);
	Create_Branches_ChargedHypotheses(locTree, locIsMCDataFlag);

	//create branches for combos
	Create_Branches_Combo(locTree, locReaction, locIsMCDataFlag, locParticleNumberMap);

	//Custom branches
	Create_CustomBranches_DataTree(locTree, locReaction, locIsMCDataFlag);
}

void DEventWriterROOT::Create_UserInfoMaps(TTree* locTree, const DReaction* locReaction, map<Particle_t, unsigned int>& locParticleNumberMap) const
{
	//kinfit type
	DKinFitType locKinFitType = locReaction->Get_KinFitType();

	//create & add reaction identification maps
	TList* locUserInfo = locTree->GetUserInfo();
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

	ostringstream locKinFitTypeStream;
	locKinFitTypeStream << locKinFitType;
	locMiscInfoMap->Add(new TObjString("KinFitType"), new TObjString(locKinFitTypeStream.str().c_str()));

	map<Particle_t, unsigned int> locParticleNumberMap_Current;
	Particle_t locPID;
	TObjString *locObjString_PID, *locObjString_Position, *locObjString_ParticleName;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		{
			ostringstream locPIDStream, locPositionStream, locParticleNameStream;
			locPID = locReactionStep->Get_InitialParticleID();
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
		}

		//target particle
		locPID = locReactionStep->Get_TargetParticleID();
		if((loc_i == 0) && (locPID != Unknown))
		{
			ostringstream locPIDStream, locPositionStream, locParticleNameStream;
			locPIDStream << PDGtype(locPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());
			locPositionStream << loc_i << "_" << -2;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());
			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
			locParticleNameStream << "Target" << Convert_ToBranchName(ParticleType(locPID));
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			string locPIDName = locParticleNameStream.str() + string("__PID");
			locMiscInfoMap->Add(new TObjString(locPIDName.c_str()), locObjString_PID);
			ostringstream locMassStream;
			locMassStream << ParticleMass(locPID);
			locMiscInfoMap->Add(new TObjString("Target__Mass"), new TObjString(locMassStream.str().c_str()));
			locParticleNameList->AddLast(locObjString_ParticleName);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			ostringstream locPIDStream, locPositionStream;
			locPID = locFinalParticleIDs[loc_j];
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

			if(locParticleNumberMap_Current.find(locPID) == locParticleNumberMap_Current.end())
				locParticleNumberMap_Current[locPID] = 1;
			else
				++locParticleNumberMap_Current[locPID];

			//see if decays in a future step
			bool locDecaysFlag = false;
			for(size_t loc_k = loc_i + 1; loc_k < locReaction->Get_NumReactionSteps(); ++loc_k)
			{
				if(locReaction->Get_ReactionStep(loc_k)->Get_InitialParticleID() != locPID)
					continue;
				locDecaysFlag = true;
				break;
			}

			ostringstream locParticleNameStream;
			if(locDecaysFlag)
				locParticleNameStream << "Decaying";
			locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
			if(locParticleNumberMap[locPID] > 1)
				locParticleNameStream << locParticleNumberMap_Current[locPID];
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
			locParticleNameList->AddLast(locObjString_ParticleName);

			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
			locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			if(locDecaysFlag)
			{
				ostringstream locMassStream;
				locMassStream << ParticleMass(locPID);
				string locMassName = locParticleNameStream.str() + string("__Mass");
				locMiscInfoMap->Add(new TObjString(locMassName.c_str()), new TObjString(locMassStream.str().c_str()));
			}
		}
	}

	//fill decay product map
	deque<size_t> locSavedSteps;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		locPID = locReactionStep->Get_InitialParticleID();
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

		//find the name of the decay parent
		//first need to find what instance this pid is a decay parent
		size_t locDecayParentInstance = 1;
		for(size_t loc_j = 0; loc_j < loc_i; ++loc_j)
		{
			if(locReaction->Get_ReactionStep(loc_j)->Get_InitialParticleID() == locPID)
				++locDecayParentInstance;
		}
		//construct the name 
		ostringstream locParticleNameStream;
		locParticleNameStream << "Decaying";
		locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
		if(locParticleNumberMap[locPID] > 1)
			locParticleNameStream << locDecayParentInstance;
		locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());

		TList* locDecayProductNames = NULL;
		Get_DecayProductNames(locReaction, loc_i, locPositionToNameMap, locDecayProductNames, locSavedSteps);
		locDecayProductMap->Add(locObjString_ParticleName, locDecayProductNames); //parent name string -> tobjarray of decay product name strings		
	}
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
		Particle_t locPID = locFinalParticleIDs[loc_j];

		//see if decays in a future step
		bool locDecaysFlag = false;
		for(size_t loc_k = locReactionStepIndex + 1; loc_k < locReaction->Get_NumReactionSteps(); ++loc_k)
		{
			if(locReaction->Get_ReactionStep(loc_k)->Get_InitialParticleID() != locPID)
				continue;
			locDecaysFlag = true;
			break;
		}

		//save and continue if doesn't decay
		if(!locDecaysFlag)
		{
			ostringstream locPositionStream;
			locPositionStream << locReactionStepIndex << "_" << loc_j;
			locDecayProductNames->AddLast(locPositionToNameMap->GetValue(locPositionStream.str().c_str()));
			continue;
		}

		//decays: find step where it decays at
		//first find what instance in the final state this particle is
		size_t locFinalStateInstance = 1;
		for(size_t loc_k = 0; loc_k < locReactionStepIndex; ++loc_k)
		{
			const DReactionStep* locNewReactionStep = locReaction->Get_ReactionStep(loc_k);
			deque<Particle_t> locNewFinalParticleIDs;
			locNewReactionStep->Get_FinalParticleIDs(locNewFinalParticleIDs);
			for(size_t loc_l = 0; loc_l < locNewFinalParticleIDs.size(); ++loc_l)
			{
				if((loc_k == locReactionStepIndex) && (loc_l == loc_j))
					break;
				if(locNewFinalParticleIDs[loc_l] == locPID)
					++locFinalStateInstance;
			}
		}

		//next find the step index that this particle decays at
		size_t locNewReactionStepIndex = 0;
		size_t locInitialStateInstance = 0;
		for(size_t loc_k = 0; loc_k < locReaction->Get_NumReactionSteps(); ++loc_k)
		{
			if(locReaction->Get_ReactionStep(loc_k)->Get_InitialParticleID() != locPID)
				continue;
			++locInitialStateInstance;
			if(locInitialStateInstance != locFinalStateInstance)
				continue;
			locNewReactionStepIndex = loc_k;
			break;
		}

		//add decay products
		Get_DecayProductNames(locReaction, locNewReactionStepIndex, locPositionToNameMap, locDecayProductNames, locSavedSteps);
	}

	locSavedSteps.push_back(locReactionStepIndex);
}

void DEventWriterROOT::Create_Branches_Thrown(TTree* locTree, bool locIsOnlyThrownFlag) const
{
	//BEAM
	Create_Branch_Fundamental<Int_t>(locTree, "ThrownBeam", "PID");
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "X4"); //reported at target center
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "P4");

	//EVENT-WIDE INFO
	Create_Branch_Fundamental<ULong64_t>(locTree, "NumPIDThrown_FinalState"); //19 digits
	Create_Branch_Fundamental<ULong64_t>(locTree, "PIDThrown_Decaying");
	Create_Branch_Fundamental<Float_t>(locTree, "MCWeight");

	//PRODUCTS
	Create_Branches_ThrownParticles(locTree, locIsOnlyThrownFlag);
}

void DEventWriterROOT::Create_Branches_ThrownParticles(TTree* locTree, bool locIsOnlyThrownFlag) const
{
	string locParticleBranchName = "Thrown";

	string locArraySizeString = "NumThrown";
	Create_Branch_Fundamental<UInt_t>(locTree, locArraySizeString);

	//IDENTIFIERS
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ParentIndex", locArraySizeString, dInitNumThrownArraySize);
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "PID", locArraySizeString, dInitNumThrownArraySize);
	if(!locIsOnlyThrownFlag)
	{
		Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "MatchID", locArraySizeString, dInitNumThrownArraySize);
		Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "MatchFOM", locArraySizeString, dInitNumThrownArraySize);
	}

	//KINEMATICS: THROWN //at the production vertex
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4", "TLorentzVector", dInitNumThrownArraySize);
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4", "TLorentzVector", dInitNumThrownArraySize);
}

void DEventWriterROOT::Create_Branches_Beam(TTree* locTree, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumBeam";
	Create_Branch_Fundamental<UInt_t>(locTree, locArraySizeString);

	//IDENTIFIER
	Create_Branch_FundamentalArray<Int_t>(locTree, "Beam", "PID", locArraySizeString, dInitNumBeamArraySize);
	if(locIsMCDataFlag)
		Create_Branch_FundamentalArray<Bool_t>(locTree, "Beam", "IsGenerator", locArraySizeString, dInitNumBeamArraySize);

	//KINEMATICS: MEASURED //at the target center
	Create_Branch_ClonesArray(locTree, "Beam", "X4_TargetCenter", "TLorentzVector", dInitNumBeamArraySize);
	Create_Branch_ClonesArray(locTree, "Beam", "P4_Measured", "TLorentzVector", dInitNumBeamArraySize);
}

void DEventWriterROOT::Create_Branches_ChargedHypotheses(TTree* locTree, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumChargedHypos";
	Create_Branch_Fundamental<UInt_t>(locTree, locArraySizeString);

	string locParticleBranchName = "ChargedHypo";

	//IDENTIFIERS / MATCHING
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "TrackID", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "PID", locArraySizeString, dInitNumTrackArraySize);
	if(locIsMCDataFlag)
		Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ThrownIndex", locArraySizeString, dInitNumTrackArraySize);

	//KINEMATICS //at the production vertex
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Measured", "TLorentzVector", dInitNumTrackArraySize);
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_Measured", "TLorentzVector", dInitNumTrackArraySize);

	//TRACKING INFO
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_Tracking", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_Tracking", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_DCdEdx", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_DCdEdx", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "dEdx_CDC", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "dEdx_FDC", locArraySizeString, dInitNumTrackArraySize);

	//TIMING INFO
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "HitTime", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "RFDeltaTVar", locArraySizeString, dInitNumTrackArraySize);

	//HIT ENERGY
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "dEdx_TOF", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "dEdx_ST", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Energy_BCAL", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Energy_FCAL", locArraySizeString, dInitNumTrackArraySize);

	//SHOWER MATCHING:
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaPhi", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaZ", locArraySizeString, dInitNumTrackArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackFCAL_DOCA", locArraySizeString, dInitNumTrackArraySize);
}

void DEventWriterROOT::Create_Branches_NeutralShowers(TTree* locTree, bool locIsMCDataFlag) const
{
	string locArraySizeString = "NumNeutralShowers";
	string locParticleBranchName = "NeutralShower";
	Create_Branch_Fundamental<UInt_t>(locTree, locArraySizeString);

	//MATCHING
	if(locIsMCDataFlag)
		Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ThrownIndex", locArraySizeString, dInitNumShowerArraySize);

	//SHOWER INFO
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Shower", "TLorentzVector", dInitNumShowerArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Energy_BCAL", locArraySizeString, dInitNumShowerArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Energy_FCAL", locArraySizeString, dInitNumShowerArraySize);

	//NEARBY TRACKS
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaPhi", locArraySizeString, dInitNumShowerArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaZ", locArraySizeString, dInitNumShowerArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "TrackFCAL_DOCA", locArraySizeString, dInitNumShowerArraySize);

	//PHOTON PID INFO
		//Computed using DVertex (best estimate of reaction vertex using all "good" tracks)
		//Can be used to compute timing chisq //is invalid for non-photons, so computed assuming photon PID
		//Variance of X4_Measured.T() - RFTime, regardless of which RF bunch is chosen. //RF bunch is combo-depende
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "PhotonRFDeltaTVar", locArraySizeString, dInitNumShowerArraySize);
}

void DEventWriterROOT::Create_Branches_Combo(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag, const map<Particle_t, unsigned int>& locParticleNumberMap) const
{
	string locNumComboString = "NumCombos";
	Create_Branch_Fundamental<UInt_t>(locTree, locNumComboString);

	//kinfit type
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);

	//create combo-dependent, particle-independent branches
	if(locIsMCDataFlag)
	{
		Create_Branch_FundamentalArray<Bool_t>(locTree, "IsTrueCombo", locNumComboString, dInitNumComboArraySize);
		Create_Branch_FundamentalArray<Bool_t>(locTree, "IsBDTSignalCombo", locNumComboString, dInitNumComboArraySize);
	}

	Create_Branch_FundamentalArray<Float_t>(locTree, "RFTime_Measured", locNumComboString, dInitNumComboArraySize);
	if(locKinFitFlag)
	{
		Create_Branch_FundamentalArray<Float_t>(locTree, "ChiSq_KinFit", locNumComboString, dInitNumComboArraySize);
		Create_Branch_FundamentalArray<UInt_t>(locTree, "NDF_KinFit", locNumComboString, dInitNumComboArraySize);
		if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
			Create_Branch_FundamentalArray<Float_t>(locTree, "RFTime_KinFit", locNumComboString, dInitNumComboArraySize);
	}

	map<Particle_t, unsigned int> locParticleNumberMap_Current;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		Particle_t locInitialPID = locReactionStep->Get_InitialParticleID();
		//should check to make sure the beam particle isn't missing...
		if((locInitialPID == Gamma) || (locInitialPID == Electron) || (locInitialPID == Positron))
			Create_Branches_BeamComboParticle(locTree, locInitialPID, locKinFitType);
		else //decaying
		{
			ostringstream locParticleNameStream;
			locParticleNameStream << "Decaying" << Convert_ToBranchName(ParticleType(locInitialPID));
			if(IsFixedMass(locInitialPID) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
				Create_Branch_ClonesArray(locTree, locParticleNameStream.str(), "P4_KinFit", "TLorentzVector", dInitNumComboArraySize);
			if(IsDetachedVertex(locInitialPID))
				Create_Branch_ClonesArray(locTree, locParticleNameStream.str(), "X4", "TLorentzVector", dInitNumComboArraySize);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			Particle_t locPID = locFinalParticleIDs[loc_j];
			//decaying particle
			if(locReaction->Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//missing particle
			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j))
			{
				// missing particle
				ostringstream locParticleNameStream;
				locParticleNameStream << "Missing" << Convert_ToBranchName(ParticleType(locPID));
				if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
					Create_Branch_ClonesArray(locTree, locParticleNameStream.str(), "P4_KinFit", "TLorentzVector", dInitNumComboArraySize);
				continue;
			}

			//detected
			if(locParticleNumberMap_Current.find(locPID) == locParticleNumberMap_Current.end())
				locParticleNumberMap_Current[locPID] = 1;
			else
				++locParticleNumberMap_Current[locPID];

			ostringstream locParticleNameStream;
			locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
			if(locParticleNumberMap.find(locPID)->second > 1)
				locParticleNameStream << locParticleNumberMap_Current[locPID];

			if(ParticleCharge(locPID) == 0)
				Create_Branches_ComboNeutral(locTree, locParticleNameStream.str(), locKinFitType);
			else
				Create_Branches_ComboTrack(locTree, locParticleNameStream.str(), locKinFitType);
		}
	}
}

void DEventWriterROOT::Create_Branches_BeamComboParticle(TTree* locTree, Particle_t locBeamPID, DKinFitType locKinFitType) const
{
	string locParticleBranchName = "ComboBeam";
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "BeamIndex", locArraySizeString, dInitNumComboArraySize);

	//SPACETIME AT INTERACTION VERTEX
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Measured", "TLorentzVector", dInitNumComboArraySize);

	//KINEMATICS: KINFIT //at the interaction vertex
	if(locKinFitType != d_NoFit)
	{
		if(((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit)) || (ParticleCharge(locBeamPID) != 0))
			Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_KinFit", "TLorentzVector", dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
			Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_KinFit", "TLorentzVector", dInitNumComboArraySize);
	}
}

void DEventWriterROOT::Create_Branches_ComboTrack(TTree* locTree, string locParticleBranchName, DKinFitType locKinFitType) const
{
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ChargedIndex", locArraySizeString, dInitNumComboArraySize); // is index to neutral/charged hypo array if PID is q0/q+-

	//PID QUALITY
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Beta_Timing_Measured", locArraySizeString, dInitNumComboArraySize);
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locArraySizeString, dInitNumComboArraySize);
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locArraySizeString, dInitNumComboArraySize);

	//KINFIT INFO //at the production vertex
	if(locKinFitType != d_NoFit)
	{
		Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_KinFit", "TLorentzVector", dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
		{
			Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_KinFit", "TLorentzVector", dInitNumComboArraySize);
			Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Beta_Timing_KinFit", locArraySizeString, dInitNumComboArraySize);
			Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locArraySizeString, dInitNumComboArraySize);
		}
	}
}

void DEventWriterROOT::Create_Branches_ComboNeutral(TTree* locTree, string locParticleBranchName, DKinFitType locKinFitType) const
{
	string locArraySizeString = "NumCombos";

	//IDENTIFIER
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ShowerIndex", locArraySizeString, dInitNumComboArraySize); // is index to neutral/charged hypo array if PID is q0/q+-

	//KINEMATICS //is combo-dependent: P4 defined by X4, X4 defined by other tracks
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Measured", "TLorentzVector", dInitNumComboArraySize);
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_Measured", "TLorentzVector", dInitNumComboArraySize);

	//PID QUALITY
	Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Beta_Timing_Measured", locArraySizeString, dInitNumComboArraySize);
	if(locParticleBranchName.substr(0, 6) == "Photon")
	{
		Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locArraySizeString, dInitNumComboArraySize);
		Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locArraySizeString, dInitNumComboArraySize);
	}

	//KINFIT INFO //at the production vertex
	if(locKinFitType != d_NoFit)
	{
		Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_KinFit", "TLorentzVector", dInitNumComboArraySize);
		if(locKinFitType != d_P4Fit)
		{
			Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_KinFit", "TLorentzVector", dInitNumComboArraySize);
			Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "Beta_Timing_KinFit", locArraySizeString, dInitNumComboArraySize);
			if(locParticleBranchName.substr(0, 6) == "Photon")
				Create_Branch_FundamentalArray<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locArraySizeString, dInitNumComboArraySize);
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

	//Pre-compute before entering lock
	ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
	Compute_ThrownPIDInfo(locMCThrowns_FinalState, locMCThrowns_Decaying, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

	//Pre-compute before entering lock
	vector<const DMCThrown*> locMCThrownsToSave;
	map<const DMCThrown*, unsigned int> locThrownIndexMap;
	Group_ThrownParticles(locMCThrowns_FinalState, locMCThrowns_Decaying, locMCThrownsToSave, locThrownIndexMap);

	string locTreeName = "Thrown_Tree";
	japp->RootWriteLock();
	{
		if(Get_ThrownTreeFileName() == "")
		{
			japp->RootUnLock();
			return;
		}

		//get the tree
		TFile* locFile = (TFile*)gROOT->FindObject(Get_ThrownTreeFileName().c_str());
		locFile->cd("/");
		TTree* locTree = (TTree*)gDirectory->Get(locTreeName.c_str());
		if(locTree == NULL)
		{
			cout << "ERROR: OUTPUT ROOT TREE NOT CREATED (in DEventWriterROOT::Fill_ThrownTree()). SKIP FILLING. " << endl;
			japp->RootUnLock();
			return;
		}

		//clear the tclonesarry's
		map<string, TClonesArray*>& locTreeMap_ClonesArray = Get_ClonesArrayMap()[locTree];
		map<string, TClonesArray*>::iterator locClonesArrayMapIterator;
		for(locClonesArrayMapIterator = locTreeMap_ClonesArray.begin(); locClonesArrayMapIterator != locTreeMap_ClonesArray.end(); ++locClonesArrayMapIterator)
			locClonesArrayMapIterator->second->Clear();

		//primary event info
		Fill_FundamentalData<UInt_t>(locTree, "RunNumber", locEventLoop->GetJEvent().GetRunNumber());
		Fill_FundamentalData<ULong64_t>(locTree, "EventNumber", locEventLoop->GetJEvent().GetEventNumber());

		//throwns
		Fill_ThrownInfo(locTree, locMCReaction, locMCThrownsToSave, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

		//Custom Branches
		Fill_CustomBranches_ThrownTree(locTree, locMCReaction, locMCThrownsToSave);

		locTree->Fill();
	}
	japp->RootUnLock();
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

	//Pre-compute before entering lock
	ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
	Compute_ThrownPIDInfo(locMCThrowns_FinalState, locMCThrowns_Decaying, locNumPIDThrown_FinalState, locPIDThrown_Decaying);

	//Pre-compute before entering lock
	vector<const DMCThrown*> locMCThrownsToSave;
	map<const DMCThrown*, unsigned int> locThrownIndexMap;
	Group_ThrownParticles(locMCThrowns_FinalState, locMCThrowns_Decaying, locMCThrownsToSave, locThrownIndexMap);

	//GET DETECTOR MATCHES
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//GET NEUTRAL INFO
	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	//GET CHARGED TRACK INFO
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	vector<const DChargedTrackHypothesis*> locComboChargedTrackHypotheses;
	locEventLoop->Get(locComboChargedTrackHypotheses, "Combo");

	//Search for hypotheses not in the PreSelect factory (i.e. PID not in REST file)
	map<const DChargedTrack*, vector<const DChargedTrackHypothesis*> > locNewHypothesesMap;
	for(size_t loc_i = 0; loc_i < locComboChargedTrackHypotheses.size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locOrigChargedTrackHypothesis = NULL;
		locComboChargedTrackHypotheses[loc_i]->GetSingle(locOrigChargedTrackHypothesis);
		if(locOrigChargedTrackHypothesis != NULL)
			continue; //PID in REST file

		const DChargedTrack* locOrigChargedTrack = NULL;
		locComboChargedTrackHypotheses[loc_i]->GetSingle(locOrigChargedTrack);

		locNewHypothesesMap[locOrigChargedTrack].push_back(locComboChargedTrackHypotheses[loc_i]);
	}

	//Build vector of combo-independent hypotheses to save
	vector<const DChargedTrackHypothesis*> locIndependentChargedTrackHypotheses;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const vector<const DChargedTrackHypothesis*>& locTrackHypoVector = locChargedTracks[loc_i]->dChargedTrackHypotheses;
		vector<const DChargedTrackHypothesis*>& locNewTrackHypoVector = locNewHypothesesMap[locChargedTracks[loc_i]];

		locIndependentChargedTrackHypotheses.insert(locIndependentChargedTrackHypotheses.end(), locTrackHypoVector.begin(), locTrackHypoVector.end());
		locIndependentChargedTrackHypotheses.insert(locIndependentChargedTrackHypotheses.end(), locNewTrackHypoVector.begin(), locNewTrackHypoVector.end());
	}

	//GET BEAM PHOTONS
		//however, only fill with beam particles that are in the combos
	set<const DBeamPhoton*> locBeamPhotonSet;
	vector<const DBeamPhoton*> locBeamPhotons;
	for(size_t loc_j = 0; loc_j < locParticleCombos.size(); ++loc_j)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombos[loc_j]->Get_ParticleComboStep(0);
		const DKinematicData* locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locKinematicData);
		if(locBeamPhoton == NULL)
			continue;
		if(locBeamPhotonSet.find(locBeamPhoton) != locBeamPhotonSet.end())
			continue; //already included
		locBeamPhotonSet.insert(locBeamPhoton);
		locBeamPhotons.push_back(locBeamPhoton);
	}

	//create map of particles to array index
		//used for pointing combo particles to the appropriate array index
	map<string, map<oid_t, int> > locObjectToArrayIndexMap;
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		locObjectToArrayIndexMap["DBeamPhoton"][locBeamPhotons[loc_i]->id] = loc_i;
	for(size_t loc_i = 0; loc_i < locIndependentChargedTrackHypotheses.size(); ++loc_i)
		locObjectToArrayIndexMap["DChargedTrackHypothesis"][locIndependentChargedTrackHypotheses[loc_i]->id] = loc_i;
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
		locObjectToArrayIndexMap["DNeutralShower"][locNeutralParticles[loc_i]->dNeutralShower->dShowerID] = loc_i;

	//EXECUTE ANALYSIS ACTIONS (Outside of ROOT Lock!)
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

	string locOutputFileName = locReaction->Get_TTreeOutputFileName();
	string locTreeName = locReaction->Get_ReactionName() + string("_Tree");

	japp->RootWriteLock();
	{
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile == NULL)
		{
			cout << "ERROR: OUTPUT ROOT TREE FILE NOT CREATED (in DEventWriterROOT::Fill_DataTree()). SKIP FILLING. " << endl;
			japp->RootUnLock();
			return;
		}

		//get the tree info
		locFile->cd("/");
		TTree* locTree = (TTree*)gDirectory->Get(locTreeName.c_str());
		if(locTree == NULL)
		{
			cout << "ERROR: OUTPUT ROOT TREE NOT CREATED (in DEventWriterROOT::Fill_DataTree()). SKIP FILLING. " << endl;
			japp->RootUnLock();
			return;
		}

		//clear the tclonesarry's
		map<string, TClonesArray*>& locTreeMap_ClonesArray = Get_ClonesArrayMap()[locTree];
		map<string, TClonesArray*>::iterator locClonesArrayMapIterator;
		for(locClonesArrayMapIterator = locTreeMap_ClonesArray.begin(); locClonesArrayMapIterator != locTreeMap_ClonesArray.end(); ++locClonesArrayMapIterator)
			locClonesArrayMapIterator->second->Clear();

		//PRIMARY EVENT INFO
		Fill_FundamentalData<UInt_t>(locTree, "RunNumber", locEventLoop->GetJEvent().GetRunNumber());
		Fill_FundamentalData<ULong64_t>(locTree, "EventNumber", locEventLoop->GetJEvent().GetEventNumber());

		//THROWN INFORMATION
		if(locMCReaction != NULL)
		{
			Fill_ThrownInfo(locTree, locMCReaction, locMCThrownsToSave, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying, locMCThrownMatching, locObjectToArrayIndexMap);
			Fill_FundamentalData<Bool_t>(locTree, "IsThrownTopology", locIsThrownTopologyFlag);
		}

		//INDEPENDENT BEAM PARTICLES
		Fill_FundamentalData<UInt_t>(locTree, "NumBeam", locBeamPhotons.size());
		for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
			Fill_BeamData(locTree, loc_i, locBeamPhotons[loc_i], locMCThrownMatching);
			//however, only fill with beam particles that are in the combos

		//INDEPENDENT CHARGED TRACKS
		Fill_FundamentalData<UInt_t>(locTree, "NumChargedHypos", locIndependentChargedTrackHypotheses.size());
		for(size_t loc_i = 0; loc_i < locIndependentChargedTrackHypotheses.size(); ++loc_i)
			Fill_ChargedHypo(locTree, loc_i, locIndependentChargedTrackHypotheses[loc_i], locMCThrownMatching, locThrownIndexMap, locDetectorMatches);

		//INDEPENDENT NEUTRAL SHOWERS
		Fill_FundamentalData<UInt_t>(locTree, "NumNeutralShowers", locNeutralParticles.size());
		for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
		{
			const DNeutralParticleHypothesis* locPhotonHypothesis = locNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
			Fill_NeutralShower(locTree, loc_i, locPhotonHypothesis, locMCThrownMatching, locThrownIndexMap, locDetectorMatches);
		}

		//COMBOS
		Fill_FundamentalData<UInt_t>(locTree, "NumCombos", locParticleCombos.size());
		for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
		{
			Fill_ComboData(locTree, locParticleCombos[loc_i], loc_i, locObjectToArrayIndexMap);
			if(locMCReaction != NULL)
			{
				Fill_FundamentalData<Bool_t>(locTree, "IsTrueCombo", locIsTrueComboFlags[loc_i], loc_i);
				Fill_FundamentalData<Bool_t>(locTree, "IsBDTSignalCombo", locIsBDTSignalComboFlags[loc_i], loc_i);
			}
		}

		//CUSTOM
		Fill_CustomBranches_DataTree(locTree, locMCReaction, locMCThrownsToSave, locMCThrownMatching, locDetectorMatches, locBeamPhotons, locIndependentChargedTrackHypotheses, locNeutralParticles, locParticleCombos);

		//FILL
		locTree->Fill();
	}
	japp->RootUnLock();
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
		unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % 10;
		if(locCurrentNumParticles != 9)
			locNumPIDThrown_FinalState += locPIDMultiplexID;
	}

	locPIDThrown_Decaying = 0;
	for(size_t loc_i = 0; loc_i < locMCThrowns_Decaying.size(); ++loc_i) //decaying
	{
		Particle_t locPID = locMCThrowns_Decaying[loc_i]->PID();
		ULong64_t locPIDMultiplexID = Calc_ParticleMultiplexID(locPID);
		if(locPID != Pi0)
			locPIDThrown_Decaying |= locPIDMultiplexID; //bit-wise or
		else //save pi0's as final state instead of decaying
		{
			unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % 10;
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

void DEventWriterROOT::Fill_ThrownInfo(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying) const
{
	map<string, map<oid_t, int> > locObjectToArrayIndexMap;
	Fill_ThrownInfo(locTree, locMCReaction, locMCThrowns, locThrownIndexMap, locNumPIDThrown_FinalState, locPIDThrown_Decaying, NULL, locObjectToArrayIndexMap);
}

void DEventWriterROOT::Fill_ThrownInfo(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying, const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	//THIS MUST BE CALLED FROM WITHIN A LOCK, SO DO NOT PASS IN JEVENTLOOP! //TOO TEMPTING TO DO SOMETHING BAD

	//WEIGHT
	Fill_FundamentalData<Float_t>(locTree, "MCWeight", locMCReaction->weight);

	//THROWN BEAM
	Fill_FundamentalData<Int_t>(locTree, "ThrownBeam", "PID", PDGtype(locMCReaction->beam.PID()));

	DVector3 locThrownBeamX3 = locMCReaction->beam.position();
	TLorentzVector locThrownBeamTX4(locThrownBeamX3.X(), locThrownBeamX3.Y(), locThrownBeamX3.Z(), locMCReaction->beam.time());
	Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "X4", locThrownBeamTX4);

	DLorentzVector locThrownBeamP4 = locMCReaction->beam.lorentzMomentum();
	TLorentzVector locThrownBeamTP4(locThrownBeamP4.Px(), locThrownBeamP4.Py(), locThrownBeamP4.Pz(), locThrownBeamP4.E());
	Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "P4", locThrownBeamTP4);

	//THROWN PRODUCTS
	Fill_FundamentalData<UInt_t>(locTree, "NumThrown", locMCThrowns.size());
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		Fill_ThrownParticleData(locTree, loc_i, locMCThrowns[loc_i], locThrownIndexMap, locMCThrownMatching, locObjectToArrayIndexMap);

	//PID INFO
	Fill_FundamentalData<ULong64_t>(locTree, "NumPIDThrown_FinalState", locNumPIDThrown_FinalState); //19 digits
	Fill_FundamentalData<ULong64_t>(locTree, "PIDThrown_Decaying", locPIDThrown_Decaying);
}

void DEventWriterROOT::Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, const DMCThrown* locMCThrown, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
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
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ParentIndex", locParentIndex, locArrayIndex);
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "PID", locMCThrown->pdgtype, locArrayIndex);

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
			//Can't use DNeutralShower JObject::id: 
				//Matching done with default-tag showers, but pre-select showers are saved to tree: JObject::id's don't match
			const DNeutralShower* locNeutralShower = locMCThrownMatching->Get_MatchingNeutralShower(locMCThrown, locMatchFOM);
			if(locNeutralShower != NULL)
				locMatchID = locObjectToArrayIndexMap.find("DNeutralShower")->second.find(locNeutralShower->dShowerID)->second;
		}
		Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "MatchID", locMatchID, locArrayIndex);
		Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "MatchFOM", locMatchFOM, locArrayIndex);
	}

	//KINEMATICS: THROWN //at the production vertex
	TLorentzVector locX4_Thrown(locMCThrown->position().X(), locMCThrown->position().Y(), locMCThrown->position().Z(), locMCThrown->time());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4", locX4_Thrown, locArrayIndex);
	TLorentzVector locP4_Thrown(locMCThrown->momentum().X(), locMCThrown->momentum().Y(), locMCThrown->momentum().Z(), locMCThrown->energy());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4", locP4_Thrown, locArrayIndex);
}

void DEventWriterROOT::Fill_BeamData(TTree* locTree, unsigned int locArrayIndex, const DBeamPhoton* locBeamPhoton, const DMCThrownMatching* locMCThrownMatching) const
{
	string locParticleBranchName = "Beam";

	//IDENTIFIER
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "PID", PDGtype(locBeamPhoton->PID()), locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Bool_t locIsGeneratorFlag = (locMCThrownMatching->Get_ReconMCGENBeamPhoton() == locBeamPhoton) ? kTRUE : kFALSE;
		Fill_FundamentalData<Bool_t>(locTree, locParticleBranchName, "IsGenerator", locIsGeneratorFlag, locArrayIndex);
	}

	//KINEMATICS: MEASURED
	DVector3 locPosition = locBeamPhoton->position();
	TLorentzVector locX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locBeamPhoton->time());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_TargetCenter", locX4_Measured, locArrayIndex);

	DLorentzVector locDP4 = locBeamPhoton->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", locP4_Measured, locArrayIndex);
}

void DEventWriterROOT::Fill_ChargedHypo(TTree* locTree, unsigned int locArrayIndex, const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const
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
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "TrackID", locChargedTrackHypothesis->candidateid, locArrayIndex);
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "PID", PDGtype(locChargedTrackHypothesis->PID()), locArrayIndex);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Int_t locThrownIndex = -1;
		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
		if(locMCThrown != NULL)
			locThrownIndex = locThrownIndexMap.find(locMCThrown)->second;
		Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ThrownIndex", locThrownIndex, locArrayIndex);
	}

	//KINEMATICS: MEASURED
	DVector3 locPosition = locChargedTrackHypothesis->position();
	TLorentzVector locTX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locChargedTrackHypothesis->time());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", locTX4_Measured, locArrayIndex);

	DLorentzVector locDP4 = locChargedTrackHypothesis->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", locP4_Measured, locArrayIndex);

	//TRACKING INFO
	Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Tracking", locChargedTrackHypothesis->dNDF_Track, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_Tracking", locChargedTrackHypothesis->dChiSq_Track, locArrayIndex);
	Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_DCdEdx", locChargedTrackHypothesis->dNDF_DCdEdx, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_DCdEdx", locChargedTrackHypothesis->dChiSq_DCdEdx, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "dEdx_CDC", locTrackTimeBased->ddEdx_CDC, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "dEdx_FDC", locTrackTimeBased->ddEdx_FDC, locArrayIndex);

	//HIT ENERGY
	double locTOFdEdx = (locChargedTrackHypothesis->Get_TOFHitMatchParams() != NULL) ? locChargedTrackHypothesis->Get_TOFHitMatchParams()->dEdx : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "dEdx_TOF", locTOFdEdx, locArrayIndex);
	double locSCdEdx = (locChargedTrackHypothesis->Get_SCHitMatchParams() != NULL) ? locChargedTrackHypothesis->Get_SCHitMatchParams()->dEdx : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "dEdx_ST", locSCdEdx, locArrayIndex);
	double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Energy_BCAL", locBCALEnergy, locArrayIndex);
	double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Energy_FCAL", locFCALEnergy, locArrayIndex);

	//TIMING INFO
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "HitTime", locChargedTrackHypothesis->t1(), locArrayIndex);
	double locStartTimeError = locChargedTrackHypothesis->t0_err();
	double locRFDeltaTVariance = (locChargedTrackHypothesis->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "RFDeltaTVar", locRFDeltaTVariance, locArrayIndex);

	//SHOWER MATCHING: BCAL
	double locTrackBCAL_DeltaPhi = 999.0, locTrackBCAL_DeltaZ = 999.0;
	if(locChargedTrackHypothesis->Get_BCALShowerMatchParams() != NULL)
	{
		locTrackBCAL_DeltaPhi = locChargedTrackHypothesis->Get_BCALShowerMatchParams()->dDeltaPhiToShower;
		locTrackBCAL_DeltaZ = locChargedTrackHypothesis->Get_BCALShowerMatchParams()->dDeltaZToShower;
	}
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaPhi", locTrackBCAL_DeltaPhi, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaZ", locTrackBCAL_DeltaZ, locArrayIndex);

	//SHOWER MATCHING: FCAL
	double locDOCAToShower_FCAL = 999.0;
	if(locChargedTrackHypothesis->Get_FCALShowerMatchParams() != NULL)
		locDOCAToShower_FCAL = locChargedTrackHypothesis->Get_FCALShowerMatchParams()->dDOCAToShower;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackFCAL_DOCA", locDOCAToShower_FCAL, locArrayIndex);
}

void DEventWriterROOT::Fill_NeutralShower(TTree* locTree, unsigned int locArrayIndex, const DNeutralParticleHypothesis* locPhotonHypothesis, const DMCThrownMatching* locMCThrownMatching, const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const
{
	string locParticleBranchName = "NeutralShower";
	const DNeutralShower* locNeutralShower = NULL;
	locPhotonHypothesis->GetSingle(locNeutralShower);

	//ASSOCIATED OBJECTS
	const DBCALShower* locBCALShower = NULL;
	locNeutralShower->GetSingle(locBCALShower);
	const DFCALShower* locFCALShower = NULL;
	locNeutralShower->GetSingle(locFCALShower);

	//MATCHING
	if(locMCThrownMatching != NULL)
	{
		Int_t locThrownIndex = -1;
		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locPhotonHypothesis, locMatchFOM);
		if(locMCThrown != NULL)
			locThrownIndex = locThrownIndexMap.find(locMCThrown)->second;
		Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ThrownIndex", locThrownIndex, locArrayIndex);
	}

	//SHOWER ENERGY
	double locBCALEnergy = (locNeutralShower->dDetectorSystem == SYS_BCAL) ? locNeutralShower->dEnergy : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Energy_BCAL", locBCALEnergy, locArrayIndex);
	double locFCALEnergy = (locNeutralShower->dDetectorSystem == SYS_FCAL) ? locNeutralShower->dEnergy : 0.0;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Energy_FCAL", locFCALEnergy, locArrayIndex);

	//SHOWER POSITION
	DLorentzVector locHitDX4 = locNeutralShower->dSpacetimeVertex;
	TLorentzVector locTX4_Shower(locHitDX4.X(), locHitDX4.Y(), locHitDX4.Z(), locHitDX4.T());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_Shower", locTX4_Shower, locArrayIndex);

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
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaPhi", locNearestTrackBCALDeltaPhi, locArrayIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackBCAL_DeltaZ", locNearestTrackBCALDeltaZ, locArrayIndex);

	//Track DOCA to Shower - FCAL
	double locDistanceToNearestTrack_FCAL = 999.0;
	if(locFCALShower != NULL)
	{
		if(!locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistanceToNearestTrack_FCAL))
			locDistanceToNearestTrack_FCAL = 999.0;
		if(locDistanceToNearestTrack_FCAL > 999.0)
			locDistanceToNearestTrack_FCAL = 999.0;
	}
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "TrackFCAL_DOCA", locDistanceToNearestTrack_FCAL, locArrayIndex);

	//PHOTON PID INFO
	double locStartTimeError = locPhotonHypothesis->t0_err();
	double PhotonRFDeltaTVar = (locPhotonHypothesis->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "PhotonRFDeltaTVar", PhotonRFDeltaTVar, locArrayIndex);
}

void DEventWriterROOT::Fill_ComboData(TTree* locTree, const DParticleCombo* locParticleCombo, unsigned int locComboIndex, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	//MAIN CLASSES
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

	//RF INFO
	double locRFTime = (locEventRFBunch != NULL) ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
	Fill_FundamentalData<Float_t>(locTree, "RFTime_Measured", locRFTime, locComboIndex);

	//KINFIT INFO
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);
	if(locKinFitFlag)
	{
		if(locKinFitResults != NULL)
		{
			Fill_FundamentalData<Float_t>(locTree, "ChiSq_KinFit", locKinFitResults->Get_ChiSq(), locComboIndex);
			Fill_FundamentalData<UInt_t>(locTree, "NDF_KinFit", locKinFitResults->Get_NDF(), locComboIndex);
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
			{
				double locRFTime_KinFit = -9.9E9; //NOT IMPLEMENTED YET
				Fill_FundamentalData<Float_t>(locTree, "RFTime_KinFit", locRFTime_KinFit, locComboIndex);
			}
		}
		else
		{
			Fill_FundamentalData<Float_t>(locTree, "ChiSq_KinFit", 0.0, locComboIndex);
			Fill_FundamentalData<UInt_t>(locTree, "NDF_KinFit", 0, locComboIndex);
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
				Fill_FundamentalData<Float_t>(locTree, "RFTime_KinFit", -9.9E9, locComboIndex);
		}
	}

	//STEP DATA
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
		Fill_ComboStepData(locTree, locParticleCombo, loc_i, locComboIndex, locKinFitType, locObjectToArrayIndexMap);
}

void DEventWriterROOT::Fill_ComboStepData(TTree* locTree, const DParticleCombo* locParticleCombo, unsigned int locStepIndex, unsigned int locComboIndex, DKinFitType locKinFitType, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const
{
	TList* locUserInfo = locTree->GetUserInfo();
	TMap* locPositionToNameMap = (TMap*)locUserInfo->FindObject("PositionToNameMap");

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
		Fill_ComboBeamData(locTree, locComboIndex, locMeasuredBeamPhoton, locBeamPhoton, locBeamIndex, locKinFitType);
	}
	else //decaying
	{
		string locParticleBranchName = string("Decaying") + Convert_ToBranchName(ParticleType(locInitialPID));
		if(IsDetachedVertex(locInitialPID))
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4", locStepTX4, locComboIndex);
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
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", locDecayP4, locComboIndex);
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

		//missing particle
		if(locParticleComboStep->Is_FinalParticleMissing(loc_i))
		{
			string locParticleBranchName = string("Missing") + Convert_ToBranchName(ParticleType(locPID));
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

				Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", locMissingP4, locComboIndex);
			}
			continue;
		}

		//get the branch name
		ostringstream locPositionStream;
		locPositionStream << locStepIndex << "_" << loc_i;
		TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
		string locParticleBranchName = (const char*)(locObjString->GetString());

		//fill the data
		if(ParticleCharge(locPID) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData);
			const DNeutralParticleHypothesis* locMeasuredNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData_Measured);
			const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
			int locShowerIndex = locObjectToArrayIndexMap.find("DNeutralShower")->second.find(locNeutralShower->dShowerID)->second;
			Fill_ComboNeutralData(locTree, locComboIndex, locParticleBranchName, locMeasuredNeutralHypo, locNeutralHypo, locShowerIndex, locKinFitType);
		}
		else
		{
			const DChargedTrackHypothesis* locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData);
			const DChargedTrackHypothesis* locMeasuredChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData_Measured);

			//find "ChargedIndex"
			const map<oid_t, int>& locObjectIDMap = locObjectToArrayIndexMap.find("DChargedTrackHypothesis")->second;
			vector<const DChargedTrackHypothesis*> locAssociatedChargedHypos;
			locKinematicData_Measured->Get(locAssociatedChargedHypos);
			int locChargedIndex = -1;
			for(size_t loc_j = 0; loc_j < locAssociatedChargedHypos.size(); ++loc_j)
			{
				map<oid_t, int>::const_iterator locIDMapIterator = locObjectIDMap.find(locAssociatedChargedHypos[loc_j]->id);
				if(locIDMapIterator == locObjectIDMap.end())
					continue;
				locChargedIndex = locIDMapIterator->second;
				break;
			}
			if(locChargedIndex == -1)
				locChargedIndex = locObjectIDMap.find(locMeasuredChargedHypo->id)->second;

			Fill_ComboChargedData(locTree, locComboIndex, locParticleBranchName, locMeasuredChargedHypo, locChargedHypo, locChargedIndex, locKinFitType);
		}
	}
}

void DEventWriterROOT::Fill_ComboBeamData(TTree* locTree, unsigned int locComboIndex, const DBeamPhoton* locMeasuredBeamPhoton, const DBeamPhoton* locBeamPhoton, unsigned int locBeamIndex, DKinFitType locKinFitType) const
{
	string locParticleBranchName = "ComboBeam";

	//IDENTIFIER
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "BeamIndex", locBeamIndex, locComboIndex);

	//MEASURED INTERACTION SPACETIME
	DVector3 locPosition = locMeasuredBeamPhoton->position();
	TLorentzVector locX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locMeasuredBeamPhoton->time());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", locX4_Measured, locComboIndex);

	//KINEMATICS: KINFIT
	if(locKinFitType != d_NoFit)
	{
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locBeamPhoton->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locBeamPhoton->time());
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", locX4_KinFit, locComboIndex);
		}

		//if charged, bends in b-field, update p4 when vertex changes
		if(((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit)) || (ParticleCharge(locMeasuredBeamPhoton->PID()) != 0))
		{
			DLorentzVector locDP4 = locBeamPhoton->lorentzMomentum();
			TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", locP4_KinFit, locComboIndex);
		}
	}
}

void DEventWriterROOT::Fill_ComboChargedData(TTree* locTree, unsigned int locComboIndex, string locParticleBranchName, const DChargedTrackHypothesis* locMeasuredChargedHypo, const DChargedTrackHypothesis* locChargedHypo, unsigned int locChargedIndex, DKinFitType locKinFitType) const
{
	//IDENTIFIER
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ChargedIndex", locChargedIndex, locComboIndex);

	//MEASURED PID INFO
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Beta_Timing_Measured", locMeasuredChargedHypo->measuredBeta(), locComboIndex);
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locMeasuredChargedHypo->dChiSq_Timing, locComboIndex);
	Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locMeasuredChargedHypo->dNDF_Timing, locComboIndex);

	//KINFIT
	if(locKinFitType != d_NoFit)
	{
		//KINEMATICS
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locChargedHypo->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locChargedHypo->time());
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", locX4_KinFit, locComboIndex);
		}

		//update even if vertex-only fit, because charged momentum propagated through b-field
		DLorentzVector locDP4 = locChargedHypo->lorentzMomentum();
		TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
		Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", locP4_KinFit, locComboIndex);

		//PID INFO
		if(locKinFitType != d_P4Fit)
		{
			Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Beta_Timing_KinFit", locChargedHypo->measuredBeta(), locComboIndex);
			Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locChargedHypo->dChiSq_Timing, locComboIndex);
		}
	}
}

void DEventWriterROOT::Fill_ComboNeutralData(TTree* locTree, unsigned int locComboIndex, string locParticleBranchName, const DNeutralParticleHypothesis* locMeasuredNeutralHypo, const DNeutralParticleHypothesis* locNeutralHypo, unsigned int locShowerIndex, DKinFitType locKinFitType) const
{
	//IDENTIFIER
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ShowerIndex", locShowerIndex, locComboIndex);

	//KINEMATICS: MEASURED
	DVector3 locPosition = locMeasuredNeutralHypo->position();
	TLorentzVector locX4_Measured(locPosition.X(), locPosition.Y(), locPosition.Z(), locMeasuredNeutralHypo->time());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", locX4_Measured, locComboIndex);

	DLorentzVector locDP4 = locMeasuredNeutralHypo->lorentzMomentum();
	TLorentzVector locP4_Measured(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
	Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", locP4_Measured, locComboIndex);

	//MEASURED PID INFO
	Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Beta_Timing_Measured", locMeasuredNeutralHypo->measuredBeta(), locComboIndex);
	if(locParticleBranchName.substr(0, 6) == "Photon")
	{
		Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locMeasuredNeutralHypo->dChiSq, locComboIndex);
		Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locMeasuredNeutralHypo->dNDF, locComboIndex);
	}

	//KINFIT
	if(locKinFitType != d_NoFit)
	{
		//KINEMATICS
		if(locKinFitType != d_P4Fit)
		{
			DVector3 locPosition = locNeutralHypo->position();
			TLorentzVector locX4_KinFit(locPosition.X(), locPosition.Y(), locPosition.Z(), locNeutralHypo->time());
			Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", locX4_KinFit, locComboIndex);
		}

		//update even if vertex-only fit, because neutral momentum defined by vertex
		DLorentzVector locDP4 = locNeutralHypo->lorentzMomentum();
		TLorentzVector locP4_KinFit(locDP4.Px(), locDP4.Py(), locDP4.Pz(), locDP4.E());
		Fill_ClonesData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", locP4_KinFit, locComboIndex);

		//PID INFO
		if(locKinFitType != d_P4Fit)
		{
			Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "Beta_Timing_KinFit", locNeutralHypo->measuredBeta(), locComboIndex);
			if(locParticleBranchName.substr(0, 6) == "Photon")
				Fill_FundamentalData<Float_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locNeutralHypo->dChiSq, locComboIndex);
		}
	}
}

