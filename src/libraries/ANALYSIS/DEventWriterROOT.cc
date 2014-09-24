#include "DEventWriterROOT.h"

unsigned int dInitNumThrownArraySize = 100;
unsigned int dInitNumUnusedArraySize = 100;
unsigned int gNumEventWriterThreads = 0;
string* dThrownTreeFileName = NULL;

map<TTree*, unsigned int>* dNumThrownArraySizeMap = NULL;
map<TTree*, unsigned int>* dNumUnusedArraySizeMap = NULL;
map<string, map<string, TClonesArray*> >* dClonesArrayMap = NULL; //first key is tree name, 2nd key is branch name
map<string, map<string, TObject*> >* dTObjectMap = NULL; //first key is tree name, 2nd key is branch name
deque<TFile*>* dOutputROOTFiles = NULL;

DEventWriterROOT::DEventWriterROOT(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	japp->RootWriteLock();
	{
		++gNumEventWriterThreads;
		if(dNumThrownArraySizeMap != NULL)
		{
			//objects already created
			japp->RootUnLock();
			return;
		}
		dNumUnusedArraySizeMap = new map<TTree*, unsigned int>();
		dNumThrownArraySizeMap = new map<TTree*, unsigned int>();
		dOutputROOTFiles = new deque<TFile*>();
		dClonesArrayMap = new map<string, map<string, TClonesArray*> >();
		dTObjectMap = new map<string, map<string, TObject*> >();
	}
	japp->RootUnLock();
}

DEventWriterROOT::~DEventWriterROOT(void)
{
	japp->RootWriteLock();
	{
		--gNumEventWriterThreads;
		if(gNumEventWriterThreads == 0)
		{
			if(dNumUnusedArraySizeMap != NULL)
			{
				delete dNumUnusedArraySizeMap;
				dNumUnusedArraySizeMap = NULL;
			}
			if(dNumThrownArraySizeMap != NULL)
			{
				delete dNumThrownArraySizeMap;
				dNumThrownArraySizeMap = NULL;
			}
			if(dOutputROOTFiles != NULL)
			{
				for(size_t loc_i = 0; loc_i < dOutputROOTFiles->size(); ++loc_i)
				{
					(*dOutputROOTFiles)[loc_i]->Write();
					(*dOutputROOTFiles)[loc_i]->Close();
					delete (*dOutputROOTFiles)[loc_i];
				}
				delete dOutputROOTFiles;
				dOutputROOTFiles = NULL;
			}
			if(dClonesArrayMap != NULL)
			{
				delete dClonesArrayMap;
				dClonesArrayMap = NULL;
			}
			if(dTObjectMap != NULL)
			{
				delete dTObjectMap;
				dTObjectMap = NULL;
			}
			if(dThrownTreeFileName != NULL)
			{
				delete dThrownTreeFileName;
				dThrownTreeFileName = NULL;
			}
		}
	}
	japp->RootUnLock();
}

void DEventWriterROOT::Create_ThrownTree(string locOutputFileName) const
{
	japp->RootWriteLock();
	{
		if(dThrownTreeFileName != NULL) //tree already created
		{
			japp->RootUnLock();
			return; //already created by another thread
		}
		dThrownTreeFileName = new string(locOutputFileName);

		//create ROOT file if it doesn't exist already
		TFile* locFile = (TFile*)gROOT->FindObject(dThrownTreeFileName->c_str());
		if(locFile == NULL)
		{
			locFile = new TFile(dThrownTreeFileName->c_str(), "RECREATE");
			dOutputROOTFiles->push_back(locFile);
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

		dNumUnusedArraySizeMap->insert(pair<TTree*, unsigned int>(locTree, 0));
		dNumThrownArraySizeMap->insert(pair<TTree*, unsigned int>(locTree, dInitNumThrownArraySize));

	/******************************************************************** Create Branches ********************************************************************/

		//create basic/misc. tree branches (run#, event#, etc.)
		Create_Branch_Fundamental<UInt_t>(locTree, "", "RunNumber", "i");
		Create_Branch_Fundamental<UInt_t>(locTree, "", "EventNumber", "i");
		Create_Branch_Fundamental<Double_t>(locTree, "", "RFTime_Thrown", "D");

		//create thrown particle branches
		//BEAM
		Create_Branch_Fundamental<Int_t>(locTree, "ThrownBeam", "PID", "I");
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "X4_Thrown", (*dTObjectMap)[locTree->GetName()]);
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "P4_Thrown", (*dTObjectMap)[locTree->GetName()]);

		//PRODUCTS
		string locNumThrownString = "NumThrown";
		Create_Branch_Fundamental<ULong64_t>(locTree, "", "NumPIDThrown_FinalState", "l"); //19 digits
		Create_Branch_Fundamental<ULong64_t>(locTree, "", "PIDThrown_Decaying", "l");
		Create_Branch_Fundamental<Double_t>(locTree, "", "MCWeight", "D");
		Create_Branch_Fundamental<UInt_t>(locTree, "", locNumThrownString, "i");
		Create_Branches_ThrownParticle(locTree, "Thrown", locNumThrownString, true);
	}
	japp->RootUnLock();
}

void DEventWriterROOT::Create_DataTrees(JEventLoop* locEventLoop) const
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

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
		dOutputROOTFiles->push_back(locFile);
	}

	//create tree if it doesn't exist already
	string locTreeName = locReactionName + string("_Tree");
	locFile->cd();
	if(gDirectory->Get(locTreeName.c_str()) != NULL)
		return; //already created by another thread
	TTree* locTree = new TTree(locTreeName.c_str(), locTreeName.c_str());

	dNumUnusedArraySizeMap->insert(pair<TTree*, unsigned int>(locTree, dInitNumUnusedArraySize));
	dNumThrownArraySizeMap->insert(pair<TTree*, unsigned int>(locTree, dInitNumThrownArraySize));

	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);

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
	map<Particle_t, unsigned int> locParticleNumberMap_Current;
	Particle_t locPID;
	ostringstream locPIDStream, locPositionStream, locParticleNameStream, locMassStream;
	TObjString *locObjString_PID, *locObjString_Position, *locObjString_ParticleName;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		locPID = locReactionStep->Get_InitialParticleID();
		locPIDStream.str("");
		locPIDStream << PDGtype(locPID);
		locObjString_PID = new TObjString(locPIDStream.str().c_str());
		locPositionStream.str("");
		locPositionStream << loc_i << "_" << -1;
		locObjString_Position = new TObjString(locPositionStream.str().c_str());
		locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
		if((loc_i == 0) && ((locPID == Gamma) || (locPID == Electron) || (locPID == Positron)))
		{
			locParticleNameStream.str("");
			locParticleNameStream << "Beam" << Convert_ToBranchName(ParticleType(locPID));
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			locParticleNameList->AddLast(locObjString_ParticleName);
		}

		//target particle
		locPID = locReactionStep->Get_TargetParticleID();
		if((loc_i == 0) && (locPID != Unknown))
		{
			locPIDStream.str("");
			locPIDStream << PDGtype(locPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());
			locPositionStream.str("");
			locPositionStream << loc_i << "_" << -2;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());
			locPositionToPIDMap->Add(locObjString_Position, locObjString_PID);
			locParticleNameStream.str("");
			locParticleNameStream << "Target" << Convert_ToBranchName(ParticleType(locPID));
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
			locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
			locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
			string locPIDName = locParticleNameStream.str() + string("__PID");
			locMiscInfoMap->Add(new TObjString(locPIDName.c_str()), locObjString_PID);
			locMassStream.str("");
			locMassStream << ParticleMass(locPID);
			locMiscInfoMap->Add(new TObjString("Target__Mass"), new TObjString(locMassStream.str().c_str()));
			locParticleNameList->AddLast(locObjString_ParticleName);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			locPID = locFinalParticleIDs[loc_j];
			locPIDStream.str("");
			locPIDStream << PDGtype(locPID);
			locObjString_PID = new TObjString(locPIDStream.str().c_str());

			locPositionStream.str("");
			locPositionStream << loc_i << "_" << loc_j;
			locObjString_Position = new TObjString(locPositionStream.str().c_str());

			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j)) //missing particle
			{
				locParticleNameStream.str("");
				locParticleNameStream << "Missing" << Convert_ToBranchName(ParticleType(locPID));
				locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());
				locNameToPositionMap->Add(locObjString_ParticleName, locObjString_Position);
				locPositionToNameMap->Add(locObjString_Position, locObjString_ParticleName);
				locNameToPIDMap->Add(locObjString_ParticleName, locObjString_PID);
				string locPIDName = locParticleNameStream.str() + string("__PID");
				locMiscInfoMap->Add(new TObjString(locPIDName.c_str()), locObjString_PID);
				locMassStream.str("");
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

			locParticleNameStream.str("");
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
				locMassStream.str("");
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
		locParticleNameStream.str("");
		locParticleNameStream << "Decaying";
		locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
		if(locParticleNumberMap[locPID] > 1)
			locParticleNameStream << locDecayParentInstance;
		locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());

		TList* locDecayProductNames = NULL;
		Get_DecayProductNames(locReaction, loc_i, locPositionToNameMap, locDecayProductNames, locSavedSteps);
		locDecayProductMap->Add(locObjString_ParticleName, locDecayProductNames); //parent name string -> tobjarray of decay product name strings		
	}

/******************************************************************** Create Branches ********************************************************************/

	//create basic/misc. tree branches (run#, event#, etc.)
	Create_Branch_Fundamental<UInt_t>(locTree, "", "RunNumber", "i");
	Create_Branch_Fundamental<UInt_t>(locTree, "", "EventNumber", "i");
	if(locIsMCDataFlag)
		Create_Branch_Fundamental<Double_t>(locTree, "", "RFTime_Thrown", "D");

	//create combo-dependent, particle-independent branches
	Create_Branch_Fundamental<Double_t>(locTree, "", "RFTime_Measured", "D");
	if(locKinFitFlag)
	{
		Create_Branch_Fundamental<Double_t>(locTree, "", "ChiSq_KinFit", "D");
		Create_Branch_Fundamental<UInt_t>(locTree, "", "NDF_KinFit", "i");
		Create_Branch_Fundamental<Double_t>(locTree, "", "RFTime_KinFit", "D");
	}

	//create thrown particle branches
	if(locIsMCDataFlag)
	{
		//BEAM
		Create_Branch_Fundamental<Int_t>(locTree, "ThrownBeam", "PID", "I");
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "X4_Thrown", (*dTObjectMap)[locTree->GetName()]);
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "ThrownBeam", "P4_Thrown", (*dTObjectMap)[locTree->GetName()]);

		//PRODUCTS
		string locNumThrownString = "NumThrown";
		Create_Branch_Fundamental<ULong64_t>(locTree, "", "NumPIDThrown_FinalState", "l"); //19 digits
		Create_Branch_Fundamental<ULong64_t>(locTree, "", "PIDThrown_Decaying", "l");
		Create_Branch_Fundamental<Double_t>(locTree, "", "MCWeight", "D");
		Create_Branch_Fundamental<UInt_t>(locTree, "", locNumThrownString, "i");
		Create_Branches_ThrownParticle(locTree, "Thrown", locNumThrownString, false);
	}

	//create branches for particles & production vertex
	locParticleNumberMap_Current.clear();
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//initial particle
		locPID = locReactionStep->Get_InitialParticleID();
		//should check to make sure the beam particle isn't missing...
		if((loc_i == 0) && ((locPID == Gamma) || (locPID == Electron) || (locPID == Positron)))
		{
			// production vertex
			Create_Branch_NoSplitTObject<TLorentzVector>(locTree, "Production", "X4", (*dTObjectMap)[locTree->GetName()]);

			// beam
			locParticleNameStream.str("");
			locParticleNameStream << "Beam" << Convert_ToBranchName(ParticleType(locPID));
			Create_Branches_Beam(locTree, locParticleNameStream.str(), locKinFitFlag);
		}
		else //decaying
		{
			locParticleNameStream.str("");
			locParticleNameStream << "Decaying" << Convert_ToBranchName(ParticleType(locPID));
			if(IsDetachedVertex(locPID))
				Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleNameStream.str(), "X4", (*dTObjectMap)[locTree->GetName()]);
			if(IsFixedMass(locPID) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
				Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleNameStream.str(), "P4", (*dTObjectMap)[locTree->GetName()]);
		}

		//final particles
		deque<Particle_t> locFinalParticleIDs;
		locReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
		for(size_t loc_j = 0; loc_j < locFinalParticleIDs.size(); ++loc_j)
		{
			locPID = locFinalParticleIDs[loc_j];
			if(locParticleNumberMap_Current.find(locPID) == locParticleNumberMap_Current.end())
				locParticleNumberMap_Current[locPID] = 1;
			else
				++locParticleNumberMap_Current[locPID];

			locParticleNameStream.str("");
			locParticleNameStream << Convert_ToBranchName(ParticleType(locPID));
			if(locParticleNumberMap[locPID] > 1)
				locParticleNameStream << locParticleNumberMap_Current[locPID];
			locObjString_ParticleName = new TObjString(locParticleNameStream.str().c_str());

			if(locReaction->Check_IsDecayingParticle(locFinalParticleIDs[loc_j], loc_i + 1))
				continue; //decaying particle
			if(locReactionStep->Get_MissingParticleIndex() == int(loc_j))
			{
				// missing particle
				locParticleNameStream.str("");
				locParticleNameStream << "Missing" << Convert_ToBranchName(ParticleType(locPID));
				Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleNameStream.str(), "X4", (*dTObjectMap)[locTree->GetName()]);
				if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
					Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleNameStream.str(), "P4", (*dTObjectMap)[locTree->GetName()]);
			}
			else //detected
				Create_Branches_FinalStateParticle(locTree, locParticleNameStream.str(), (ParticleCharge(locPID) != 0), locKinFitFlag, locIsMCDataFlag);
		}
	}

	//create unused particle branches
	string locNumUnusedString = "NumUnused";
	Create_Branch_Fundamental<UInt_t>(locTree, "", locNumUnusedString, "i");
	Create_Branches_UnusedParticle(locTree, "Unused", locNumUnusedString, locIsMCDataFlag);
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

void DEventWriterROOT::Create_Branches_FinalStateParticle(TTree* locTree, string locParticleBranchName, bool locIsChargedFlag, bool locKinFitFlag, bool locIsMCDataFlag) const
{
	//IDENTIFIER / MATCHING
	Create_Branch_Fundamental<Int_t>(locTree, locParticleBranchName, "ObjectID", "I");
	if(locIsMCDataFlag)
		Create_Branch_Fundamental<Int_t>(locTree, locParticleBranchName, "MatchID", "I");

	//KINEMATICS: MEASURED //at the production vertex
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", (*dTObjectMap)[locTree->GetName()]);
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", (*dTObjectMap)[locTree->GetName()]);

	//KINEMATICS: KINFIT //at the production vertex
	if(locKinFitFlag)
	{
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", (*dTObjectMap)[locTree->GetName()]);
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", (*dTObjectMap)[locTree->GetName()]);
	}

	//PID QUALITY
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "AvgBeta_Timing", "D");
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "HitTime", "D");
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", "D");
	if(locKinFitFlag)
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", "D");
	Create_Branch_Fundamental<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", "i");
	if(locIsChargedFlag)
	{
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "ChiSq_Tracking", "D");
		Create_Branch_Fundamental<UInt_t>(locTree, locParticleBranchName, "NDF_Tracking", "i");
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "ChiSq_DCdEdx", "D");
		Create_Branch_Fundamental<UInt_t>(locTree, locParticleBranchName, "NDF_DCdEdx", "i");
	}

	//DEPOSITED ENERGY
	if(locIsChargedFlag)
	{
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "dEdx_CDC", "D");
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "dEdx_FDC", "D");
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "dEdx_TOF", "D");
		Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "dEdx_ST", "D");
	}
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "Energy_BCAL", "D");
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "Energy_FCAL", "D");
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "TrackDOCA_BCAL", "D");
	Create_Branch_Fundamental<Double_t>(locTree, locParticleBranchName, "TrackDOCA_FCAL", "D");
}

void DEventWriterROOT::Create_Branches_Beam(TTree* locTree, string locParticleBranchName, bool locKinFitFlag) const
{
	//IDENTIFIER
	Create_Branch_Fundamental<Int_t>(locTree, locParticleBranchName, "ObjectID", "I");

	//KINEMATICS: MEASURED //at the target center
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", (*dTObjectMap)[locTree->GetName()]);
	Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", (*dTObjectMap)[locTree->GetName()]);

	//KINEMATICS: KINFIT //at the interaction vertex
	if(locKinFitFlag)
	{
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", (*dTObjectMap)[locTree->GetName()]);
		Create_Branch_NoSplitTObject<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", (*dTObjectMap)[locTree->GetName()]);
	}
}

void DEventWriterROOT::Create_Branches_UnusedParticle(TTree* locTree, string locParticleBranchName, string locArraySizeString, bool locIsMCDataFlag) const
{
	//IDENTIFIERS / MATCHING
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ObjectID", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "I");
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "PID", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "I");
	if(locIsMCDataFlag)
		Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "MatchID", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "I");

	//KINEMATICS: MEASURED //at the production vertex
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Measured", "TLorentzVector", (*dNumUnusedArraySizeMap)[locTree]);
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_Measured", "TLorentzVector", (*dNumUnusedArraySizeMap)[locTree]);

	//PID QUALITY
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "ChiSq_Tracking", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_Tracking", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "i");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "AvgBeta_Timing", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "HitTime", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "i");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "ChiSq_DCdEdx", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<UInt_t>(locTree, locParticleBranchName, "NDF_DCdEdx", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "i");

	//DEPOSITED ENERGY
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "dEdx_CDC", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "dEdx_FDC", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "dEdx_TOF", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "dEdx_ST", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "Energy_BCAL", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "Energy_FCAL", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "TrackDOCA_BCAL", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
	Create_Branch_FundamentalArray<Double_t>(locTree, locParticleBranchName, "TrackDOCA_FCAL", locArraySizeString, (*dNumUnusedArraySizeMap)[locTree], "D");
}

void DEventWriterROOT::Create_Branches_ThrownParticle(TTree* locTree, string locParticleBranchName, string locArraySizeString, bool locIsOnlyThrownFlag) const
{
	//IDENTIFIERS
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "ParentID", locArraySizeString, (*dNumThrownArraySizeMap)[locTree], "I");
	Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "PID", locArraySizeString, (*dNumThrownArraySizeMap)[locTree], "I");
	if(!locIsOnlyThrownFlag)
		Create_Branch_FundamentalArray<Int_t>(locTree, locParticleBranchName, "MatchID", locArraySizeString, (*dNumThrownArraySizeMap)[locTree], "I");

	//KINEMATICS: THROWN //at the production vertex
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "X4_Thrown", "TLorentzVector", (*dNumThrownArraySizeMap)[locTree]);
	Create_Branch_ClonesArray(locTree, locParticleBranchName, "P4_Thrown", "TLorentzVector", (*dNumThrownArraySizeMap)[locTree]);
}

string DEventWriterROOT::Create_Branch_ClonesArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locClassName, unsigned int locSize) const
{
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	(*dClonesArrayMap)[locTree->GetName()].insert(pair<string, TClonesArray*>(locBranchName, new TClonesArray(locClassName.c_str(), locSize)));
	locTree->Branch(locBranchName.c_str(), &((*dClonesArrayMap)[locTree->GetName()][locBranchName]), 32000, 0); //0: don't split
	return locBranchName;
}

string DEventWriterROOT::Convert_ToBranchName(string locInputName) const
{
	TString locTString(locInputName);
	locTString.ReplaceAll("*", "Star");
	locTString.ReplaceAll("(", "_");
	locTString.ReplaceAll(")", "_");
	locTString.ReplaceAll("+", "Plus");
	locTString.ReplaceAll("-", "Minus");
	locTString.ReplaceAll("'", "Prime");
	return (string)((const char*)locTString);
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
	if(dThrownTreeFileName == NULL)
		return;

	vector<const DMCThrown*> locMCThrowns_FinalState;
	locEventLoop->Get(locMCThrowns_FinalState, "FinalState");

	vector<const DMCThrown*> locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");

	vector<const DEventRFBunch*> locThrownEventRFBunches;
	locEventLoop->Get(locThrownEventRFBunches, "Thrown");

	const DMCReaction* locMCReaction = NULL;
	locEventLoop->GetSingle(locMCReaction);

	//create map of mcthrown to array index
	map<const DMCThrown*, unsigned int> locThrownObjectIDMap;
	for(size_t loc_i = 0; loc_i < locMCThrowns_FinalState.size(); ++loc_i)
		locThrownObjectIDMap[locMCThrowns_FinalState[loc_i]] = loc_i;
	for(size_t loc_i = 0; loc_i < locMCThrowns_Decaying.size(); ++loc_i)
		locThrownObjectIDMap[locMCThrowns_Decaying[loc_i]] = loc_i + locMCThrowns_FinalState.size();

	string locTreeName = "Thrown_Tree";
	japp->RootWriteLock();
	{
		//get the tree
		TFile* locFile = (TFile*)gROOT->FindObject(dThrownTreeFileName->c_str());
		locFile->cd("/");
		TTree* locTree = (TTree*)gDirectory->Get(locTreeName.c_str());
		if(locTree == NULL)
		{
			cout << "ERROR: OUTPUT ROOT TREE NOT CREATED (in DEventWriterROOT::Fill_ThrownTree()). SKIP FILLING. " << endl;
			japp->RootUnLock();
			return;
		}

		//clear the tclonesarry's
		map<string, TClonesArray*>& locTreeMap_ClonesArray = (*dClonesArrayMap)[locTreeName];
		map<string, TClonesArray*>::iterator locClonesArrayMapIterator;
		for(locClonesArrayMapIterator = locTreeMap_ClonesArray.begin(); locClonesArrayMapIterator != locTreeMap_ClonesArray.end(); ++locClonesArrayMapIterator)
			locClonesArrayMapIterator->second->Clear();

		//primary event info
		Fill_FundamentalData<UInt_t>(locTree, "RunNumber", locEventLoop->GetJEvent().GetRunNumber());
		Fill_FundamentalData<UInt_t>(locTree, "EventNumber", locEventLoop->GetJEvent().GetEventNumber());
		if(!locThrownEventRFBunches.empty())
			Fill_FundamentalData<Double_t>(locTree, "RFTime_Thrown", locThrownEventRFBunches[0]->dTime);

		//throwns
		size_t locNumThrown = locMCThrowns_FinalState.size() + locMCThrowns_Decaying.size();
		if(locNumThrown > 0)
		{
			Double_t locMCWeight = (locMCReaction == NULL) ? 0.0 : locMCReaction->weight;
			Fill_FundamentalData<Double_t>(locTree, "MCWeight", locMCWeight);

			//THROWN BEAM
			Fill_FundamentalData<Int_t>(locTree, "ThrownBeam__PID", PDGtype(locMCReaction->beam.PID()));

			DVector3 locThrownBeamX3 = locMCReaction->beam.position();
			TLorentzVector locThrownBeamTX4(locThrownBeamX3.X(), locThrownBeamX3.Y(), locThrownBeamX3.Z(), locMCReaction->beam.time());
			Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "X4_Thrown", &locThrownBeamTX4, (*dTObjectMap)[locTree->GetName()]);

			DLorentzVector locThrownBeamP4 = locMCReaction->beam.lorentzMomentum();
			TLorentzVector locThrownBeamTP4(locThrownBeamP4.Px(), locThrownBeamP4.Py(), locThrownBeamP4.Pz(), locThrownBeamP4.E());
			Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "P4_Thrown", &locThrownBeamTP4, (*dTObjectMap)[locTree->GetName()]);

			//THROWN PRODUCTS
			Fill_FundamentalData<UInt_t>(locTree, "NumThrown", locNumThrown);
			for(size_t loc_j = 0; loc_j < locMCThrowns_FinalState.size(); ++loc_j)
				Fill_ThrownParticleData(locTree, loc_j, locNumThrown, locMCThrowns_FinalState[loc_j], locThrownObjectIDMap);
			for(size_t loc_j = 0; loc_j < locMCThrowns_Decaying.size(); ++loc_j)
				Fill_ThrownParticleData(locTree, loc_j + locMCThrowns_FinalState.size(), locNumThrown, locMCThrowns_Decaying[loc_j], locThrownObjectIDMap);

			//THROWN PARTICLES BY PID
			ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
			for(size_t loc_j = 0; loc_j < locMCThrowns_FinalState.size(); ++loc_j) //final state
			{
				Particle_t locPID = locMCThrowns_FinalState[loc_j]->PID();
				ULong64_t locPIDMultiplexID = Calc_ParticleMultiplexID(locPID);
				unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % 10;
				if(locCurrentNumParticles != 9)
					locNumPIDThrown_FinalState += locPIDMultiplexID;
			}
			for(size_t loc_j = 0; loc_j < locMCThrowns_Decaying.size(); ++loc_j) //decaying
			{
				Particle_t locPID = locMCThrowns_Decaying[loc_j]->PID();
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

			Fill_FundamentalData<ULong64_t>(locTree, "NumPIDThrown_FinalState", locNumPIDThrown_FinalState); //19 digits
			Fill_FundamentalData<ULong64_t>(locTree, "PIDThrown_Decaying", locPIDThrown_Decaying);

			//update array sizes
			if(locNumThrown > (*dNumThrownArraySizeMap)[locTree])
				(*dNumThrownArraySizeMap)[locTree] = locNumThrown;
		}

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

	vector<const DMCThrown*> locMCThrowns_FinalState;
	locEventLoop->Get(locMCThrowns_FinalState, "FinalState");

	vector<const DMCThrown*> locMCThrowns_Decaying;
	locEventLoop->Get(locMCThrowns_Decaying, "Decaying");

	vector<const DEventRFBunch*> locThrownEventRFBunches;
	locEventLoop->Get(locThrownEventRFBunches, "Thrown");

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses);

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

	const DAnalysisUtilities* locAnalysisUtilities = NULL;
	locEventLoop->GetSingle(locAnalysisUtilities);

	const DMCReaction* locMCReaction = NULL;
	locEventLoop->GetSingle(locMCReaction);

	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	bool locKinFitFlag = (locKinFitType != d_NoFit);

	//find max charged identifier #:
	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);
	unsigned long locMaxChargedID = 0;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		unsigned long locCandidateID = locChargedTracks[loc_i]->Get_BestFOM()->candidateid;
		if(locCandidateID > locMaxChargedID)
			locMaxChargedID = locCandidateID;
	}

	//create map of neutral shower to new id #
	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers);
	map<const DNeutralShower*, int> locShowerToIDMap;
	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
		locShowerToIDMap[locNeutralShowers[loc_i]] = locMaxChargedID + 1 + loc_i;

	//create map of beam photon to new id #
	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons);
	map<const DBeamPhoton*, int> locBeamToIDMap;
	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); ++loc_i)
		locBeamToIDMap[locBeamPhotons[loc_i]] = loc_i;

	//create map of mcthrown to array index
	map<const DMCThrown*, unsigned int> locThrownObjectIDMap;
	for(size_t loc_i = 0; loc_i < locMCThrowns_FinalState.size(); ++loc_i)
		locThrownObjectIDMap[locMCThrowns_FinalState[loc_i]] = loc_i;
	for(size_t loc_i = 0; loc_i < locMCThrowns_Decaying.size(); ++loc_i)
		locThrownObjectIDMap[locMCThrowns_Decaying[loc_i]] = loc_i + locMCThrowns_FinalState.size();

	double locMinThrownMatchFOM = locReaction->Get_MinThrownMatchFOMForROOT();
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
		TList* locUserInfo = locTree->GetUserInfo();
		TMap* locPositionToNameMap = (TMap*)locUserInfo->FindObject("PositionToNameMap");

		//loop over combos (tree entries)
		for(size_t loc_j = 0; loc_j < locParticleCombos.size(); ++loc_j)
		{
			//clear the tclonesarry's
			map<string, TClonesArray*>& locTreeMap_ClonesArray = (*dClonesArrayMap)[locTreeName];
			map<string, TClonesArray*>::iterator locClonesArrayMapIterator;
			for(locClonesArrayMapIterator = locTreeMap_ClonesArray.begin(); locClonesArrayMapIterator != locTreeMap_ClonesArray.end(); ++locClonesArrayMapIterator)
				locClonesArrayMapIterator->second->Clear();

			//primary event info
			Fill_FundamentalData<UInt_t>(locTree, "RunNumber", locEventLoop->GetJEvent().GetRunNumber());
			Fill_FundamentalData<UInt_t>(locTree, "EventNumber", locEventLoop->GetJEvent().GetEventNumber());
			if(!locThrownEventRFBunches.empty())
				Fill_FundamentalData<Double_t>(locTree, "RFTime_Thrown", locThrownEventRFBunches[0]->dTime);

			//throwns
			size_t locNumThrown = locMCThrowns_FinalState.size() + locMCThrowns_Decaying.size();
			if(locNumThrown > 0)
			{
				Double_t locMCWeight = (locMCReaction == NULL) ? 0.0 : locMCReaction->weight;
				Fill_FundamentalData<Double_t>(locTree, "MCWeight", locMCWeight);

				//THROWN BEAM
				Fill_FundamentalData<Int_t>(locTree, "ThrownBeam__PID", PDGtype(locMCReaction->beam.PID()));

				DVector3 locThrownBeamX3 = locMCReaction->beam.position();
				TLorentzVector locThrownBeamTX4(locThrownBeamX3.X(), locThrownBeamX3.Y(), locThrownBeamX3.Z(), locMCReaction->beam.time());
				Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "X4_Thrown", &locThrownBeamTX4, (*dTObjectMap)[locTree->GetName()]);

				DLorentzVector locThrownBeamP4 = locMCReaction->beam.lorentzMomentum();
				TLorentzVector locThrownBeamTP4(locThrownBeamP4.Px(), locThrownBeamP4.Py(), locThrownBeamP4.Pz(), locThrownBeamP4.E());
				Fill_TObjectData<TLorentzVector>(locTree, "ThrownBeam", "P4_Thrown", &locThrownBeamTP4, (*dTObjectMap)[locTree->GetName()]);

				//THROWN PRODUCTS
				Fill_FundamentalData<UInt_t>(locTree, "NumThrown", locNumThrown);
				for(size_t loc_j = 0; loc_j < locMCThrowns_FinalState.size(); ++loc_j)
					Fill_ThrownParticleData(locTree, loc_j, locNumThrown, locMCThrowns_FinalState[loc_j], locThrownObjectIDMap, locMCThrownMatching, locMinThrownMatchFOM, locShowerToIDMap);
				for(size_t loc_j = 0; loc_j < locMCThrowns_Decaying.size(); ++loc_j)
					Fill_ThrownParticleData(locTree, loc_j + locMCThrowns_FinalState.size(), locNumThrown, locMCThrowns_Decaying[loc_j], locThrownObjectIDMap, locMCThrownMatching, locMinThrownMatchFOM, locShowerToIDMap);

				//THROWN PARTICLES BY PID
				ULong64_t locNumPIDThrown_FinalState = 0, locPIDThrown_Decaying = 0;
				for(size_t loc_j = 0; loc_j < locMCThrowns_FinalState.size(); ++loc_j) //final state
				{
					Particle_t locPID = locMCThrowns_FinalState[loc_j]->PID();
					ULong64_t locPIDMultiplexID = Calc_ParticleMultiplexID(locPID);
					unsigned int locCurrentNumParticles = (locNumPIDThrown_FinalState / locPIDMultiplexID) % 10;
					if(locCurrentNumParticles != 9)
						locNumPIDThrown_FinalState += locPIDMultiplexID;
				}
				for(size_t loc_j = 0; loc_j < locMCThrowns_Decaying.size(); ++loc_j) //decaying
				{
					Particle_t locPID = locMCThrowns_Decaying[loc_j]->PID();
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

				Fill_FundamentalData<ULong64_t>(locTree, "NumPIDThrown_FinalState", locNumPIDThrown_FinalState); //19 digits
				Fill_FundamentalData<ULong64_t>(locTree, "PIDThrown_Decaying", locPIDThrown_Decaying);

				//update array sizes
				if(locNumThrown > (*dNumThrownArraySizeMap)[locTree])
					(*dNumThrownArraySizeMap)[locTree] = locNumThrown;
			}

			//combo data
			const DParticleCombo* locParticleCombo = locParticleCombos[loc_j];
			const DKinFitResults* locKinFitResults = locParticleCombo->Get_KinFitResults();
			const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

			//rf & kinfit data
			double locRFTime = (locEventRFBunch != NULL) ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
			Fill_FundamentalData<Double_t>(locTree, "RFTime_Measured", locRFTime);
			if(locKinFitFlag)
			{
				if(locKinFitResults != NULL)
				{
					Fill_FundamentalData<Double_t>(locTree, "ChiSq_KinFit", locKinFitResults->Get_ChiSq());
					Fill_FundamentalData<UInt_t>(locTree, "NDF_KinFit", locKinFitResults->Get_NDF());
				}
				else
				{
					Fill_FundamentalData<Double_t>(locTree, "ChiSq_KinFit", 0.0);
					Fill_FundamentalData<UInt_t>(locTree, "NDF_KinFit", 0);
				}
				double locRFTime_KinFit = (locEventRFBunch != NULL) ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();
				Fill_FundamentalData<Double_t>(locTree, "RFTime_KinFit", locRFTime_KinFit);
			}

			//steps
			vector<const DChargedTrackHypothesis*> locUnusedChargedTrackHypotheses = locChargedTrackHypotheses;
			vector<const DNeutralParticleHypothesis*> locUnusedNeutralParticleHypotheses = locNeutralParticleHypotheses;
			for(size_t loc_k = 0; loc_k < locParticleCombo->Get_NumParticleComboSteps(); ++loc_k)
			{
				const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_k);
				DLorentzVector locStepX4 = locParticleComboStep->Get_SpacetimeVertex();
				TLorentzVector locStepTX4(locStepX4.X(), locStepX4.Y(), locStepX4.Z(), locStepX4.T());

				//beam & production vertex
				Particle_t locPID = locParticleComboStep->Get_InitialParticleID();
				const DKinematicData* locKinematicData = locParticleComboStep->Get_InitialParticle();
				if((loc_k == 0) && ((locPID == Gamma) || (locPID == Electron) || (locPID == Positron)))
				{
					const DKinematicData* locKinematicData_Measured = locParticleComboStep->Get_InitialParticle_Measured();
					string locParticleBranchName = string("Beam") + Convert_ToBranchName(ParticleType(locPID));
					if(locKinematicData_Measured != NULL) //missing beam particle
						Fill_BeamParticleData(locTree, locParticleBranchName, locKinematicData, locKinematicData_Measured, locBeamToIDMap);

					// production vertex
					Fill_TObjectData<TLorentzVector>(locTree, "Production", "X4", &locStepTX4, (*dTObjectMap)[locTree->GetName()]);
				}
				else //decaying
				{
					string locParticleBranchName = string("Decaying") + Convert_ToBranchName(ParticleType(locPID));
					if(IsDetachedVertex(locPID))
						Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4", &locStepTX4, (*dTObjectMap)[locTree->GetName()]);
					if(IsFixedMass(locPID) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
					{
						TLorentzVector locDecayP4;
						if(locKinFitResults == NULL)
						{
							//fit failed to converge, calc from other particles
							DLorentzVector locDecayDP4 = locAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, loc_k, false);
							locDecayDP4.SetE(sqrt(locDecayDP4.Vect().Mag2() + ParticleMass(locPID)*ParticleMass(locPID)));
							locDecayP4.SetPxPyPzE(locDecayDP4.Px(), locDecayDP4.Py(), locDecayDP4.Pz(), locDecayDP4.E());
						}
						else
							locDecayP4.SetPxPyPzE(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());
						Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4", &locDecayP4, (*dTObjectMap)[locTree->GetName()]);
					}
				}

				//final state particles
				for(size_t loc_l = 0; loc_l < locParticleComboStep->Get_NumFinalParticles(); ++loc_l)
				{
					locPID = locParticleComboStep->Get_FinalParticleID(loc_l);
					const DKinematicData* locKinematicData = locParticleComboStep->Get_FinalParticle(loc_l);
					const DKinematicData* locKinematicData_Measured = locParticleComboStep->Get_FinalParticle_Measured(loc_l);

					//decaying particle
					if(locParticleComboStep->Is_FinalParticleDecaying(loc_l))
						continue;

					//missing particle
					if(locParticleComboStep->Is_FinalParticleMissing(loc_l))
					{
						string locParticleBranchName = string("Missing") + Convert_ToBranchName(ParticleType(locPID));
						Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4", &locStepTX4, (*dTObjectMap)[locTree->GetName()]);
						if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
						{
							TLorentzVector locMissingP4;
							if(locKinFitResults == NULL)
							{
								//fit failed to converge, calc from other particles
								DLorentzVector locMissingDP4 = locAnalysisUtilities->Calc_MissingP4(locParticleCombo, false);
								locMissingDP4.SetE(sqrt(locMissingDP4.Vect().Mag2() + ParticleMass(locPID)*ParticleMass(locPID)));
								locMissingP4.SetPxPyPzE(locMissingDP4.Px(), locMissingDP4.Py(), locMissingDP4.Pz(), locMissingDP4.E());
							}
							else
								locMissingP4.SetPxPyPzE(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());
							Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4", &locMissingP4, (*dTObjectMap)[locTree->GetName()]);
						}
						continue;
					}

					//remove from unused lists
					const JObject* locSourceObject = locParticleComboStep->Get_FinalParticle_SourceObject(loc_l);
					if(locParticleComboStep->Is_FinalParticleCharged(loc_l))
					{
						const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locSourceObject);
						oid_t locTrackID = locChargedTrack->Get_BestFOM()->candidateid;
						for(size_t loc_m = 0; loc_m < locUnusedChargedTrackHypotheses.size(); ++loc_m)
						{
							if((locUnusedChargedTrackHypotheses[loc_m]->candidateid != locTrackID) || (locUnusedChargedTrackHypotheses[loc_m]->PID() != locKinematicData_Measured->PID()))
								continue;
							locUnusedChargedTrackHypotheses.erase(locUnusedChargedTrackHypotheses.begin() + loc_m);
							break;
						}
					}
					else
					{
						const DNeutralShower* locSourceNeutralShower = dynamic_cast<const DNeutralShower*>(locSourceObject);
						for(size_t loc_m = 0; loc_m < locUnusedNeutralParticleHypotheses.size(); ++loc_m)
						{
							const DNeutralShower* locAssociatedNeutralShower = NULL;
							locUnusedNeutralParticleHypotheses[loc_m]->GetSingleT(locAssociatedNeutralShower);
							if((locSourceNeutralShower != locAssociatedNeutralShower) || (locUnusedNeutralParticleHypotheses[loc_m]->PID() != locKinematicData_Measured->PID()))
								continue;
							locUnusedNeutralParticleHypotheses.erase(locUnusedNeutralParticleHypotheses.begin() + loc_m);
							break;
						}
					}

					//get the branch name
					ostringstream locPositionStream;
					locPositionStream << loc_k << "_" << loc_l;
					TObjString* locObjString = (TObjString*)locPositionToNameMap->GetValue(locPositionStream.str().c_str());
					string locParticleBranchName = (const char*)(locObjString->GetString());
					//fill the data
					Fill_ParticleData(locKinFitFlag, locTree, locParticleBranchName, locKinematicData, locKinematicData_Measured, locEventRFBunch, locShowerToIDMap, locMCThrownMatching, locMinThrownMatchFOM, locThrownObjectIDMap, locDetectorMatches);
				}
			}

			//unused
			unsigned int locNumUnused = locUnusedChargedTrackHypotheses.size() + locUnusedNeutralParticleHypotheses.size();
			Fill_FundamentalData<UInt_t>(locTree, "NumUnused", locNumUnused);
			for(size_t loc_j = 0; loc_j < locUnusedChargedTrackHypotheses.size(); ++loc_j)
				Fill_UnusedParticleData(locTree, loc_j, locNumUnused, locUnusedChargedTrackHypotheses[loc_j], locEventRFBunch, locShowerToIDMap, locMCThrownMatching, locMinThrownMatchFOM, locThrownObjectIDMap, locDetectorMatches);
			for(size_t loc_j = 0; loc_j < locUnusedNeutralParticleHypotheses.size(); ++loc_j)
				Fill_UnusedParticleData(locTree, loc_j + locUnusedChargedTrackHypotheses.size(), locNumUnused, locUnusedNeutralParticleHypotheses[loc_j], locEventRFBunch, locShowerToIDMap, locMCThrownMatching, locMinThrownMatchFOM, locThrownObjectIDMap, locDetectorMatches);
			//update array sizes
			if(locNumUnused > (*dNumUnusedArraySizeMap)[locTree])
				(*dNumUnusedArraySizeMap)[locTree] = locNumUnused;

			locTree->Fill();
		}
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

void DEventWriterROOT::Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DMCThrown* locMCThrown, map<const DMCThrown*, unsigned int> locThrownObjectIDMap) const
{
	map<const DNeutralShower*, int> locShowerToIDMap;
	Fill_ThrownParticleData(locTree, locArrayIndex, locMinArraySize, locMCThrown, locThrownObjectIDMap, NULL, 0.0, locShowerToIDMap);
}

void DEventWriterROOT::Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DMCThrown* locMCThrown, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, const map<const DNeutralShower*, int>& locShowerToIDMap) const
{
	//IDENTIFIERS
	int locParentID = -1; //e.g. photoproduced
	map<const DMCThrown*, unsigned int>::const_iterator locIterator;
	for(locIterator = locThrownObjectIDMap.begin(); locIterator != locThrownObjectIDMap.end(); ++locIterator)
	{
		if(locIterator->first->myid != locMCThrown->parentid)
			continue;
		locParentID = locIterator->second;
		break;
	}
	Fill_FundamentalData<Int_t>(locTree, "Thrown", "ParentID", locParentID, locArrayIndex, locMinArraySize, (*dNumThrownArraySizeMap)[locTree]);
	Fill_FundamentalData<Int_t>(locTree, "Thrown", "PID", locMCThrown->pdgtype, locArrayIndex, locMinArraySize, (*dNumThrownArraySizeMap)[locTree]);
	Int_t locMatchID = -1;
	if(locMCThrownMatching != NULL)
	{
		if(ParticleCharge(locMCThrown->PID()) != 0)
		{
			double locMatchFOM = 0.0;
			const DChargedTrack* locChargedTrack = locMCThrownMatching->Get_MatchingChargedTrack(locMCThrown, locMatchFOM);
			if((locChargedTrack != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
				locMatchID = locChargedTrack->Get_BestFOM()->candidateid;
		}
		else
		{
			double locMatchFOM = 0.0;
			const DNeutralParticle* locNeutralParticle = locMCThrownMatching->Get_MatchingNeutralParticle(locMCThrown, locMatchFOM);
			if((locNeutralParticle != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
			{
				const DNeutralShower* locNeutralShower = NULL;
				locNeutralParticle->GetSingleT(locNeutralShower);
				locMatchID = (locShowerToIDMap.find(locNeutralShower))->second;
			}
		}
		Fill_FundamentalData<Int_t>(locTree, "Thrown", "MatchID", locMatchID, locArrayIndex, locMinArraySize, (*dNumThrownArraySizeMap)[locTree]);
	}

	//KINEMATICS: THROWN //at the production vertex
	TLorentzVector locX4_Thrown(locMCThrown->position().X(), locMCThrown->position().Y(), locMCThrown->position().Z(), locMCThrown->time());
	Fill_ClonesData<TLorentzVector>(locTree, "Thrown", "X4_Thrown", &locX4_Thrown, locArrayIndex, (*dClonesArrayMap)[locTree->GetName()]);
	TLorentzVector locP4_Thrown(locMCThrown->momentum().X(), locMCThrown->momentum().Y(), locMCThrown->momentum().Z(), locMCThrown->energy());
	Fill_ClonesData<TLorentzVector>(locTree, "Thrown", "P4_Thrown", &locP4_Thrown, locArrayIndex, (*dClonesArrayMap)[locTree->GetName()]);
}

void DEventWriterROOT::Fill_BeamParticleData(TTree* locTree, string locParticleBranchName, const DKinematicData* locKinematicData, const DKinematicData* locKinematicData_Measured, const map<const DBeamPhoton*, int>& locBeamToIDMap) const
{
	//IDENTIFIER
	const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locKinematicData_Measured);
	Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ObjectID", locBeamToIDMap.find(locBeamPhoton)->second);

	//KINEMATICS: MEASURED
	TLorentzVector locX4_Measured(locKinematicData_Measured->position().X(), locKinematicData_Measured->position().Y(), locKinematicData_Measured->position().Z(), locKinematicData_Measured->time());
	Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", &locX4_Measured, (*dTObjectMap)[locTree->GetName()]);
	TLorentzVector locP4_Measured(locKinematicData_Measured->momentum().X(), locKinematicData_Measured->momentum().Y(), locKinematicData_Measured->momentum().Z(), locKinematicData_Measured->energy());
	Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", &locP4_Measured, (*dTObjectMap)[locTree->GetName()]);

	//KINEMATICS: KINFIT
	if(locKinematicData != locKinematicData_Measured)
	{
		TLorentzVector locX4_KinFit(locKinematicData->position().X(), locKinematicData->position().Y(), locKinematicData->position().Z(), locKinematicData->time());
		Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", &locX4_KinFit, (*dTObjectMap)[locTree->GetName()]);
		TLorentzVector locP4_KinFit(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());
		Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", &locP4_KinFit, (*dTObjectMap)[locTree->GetName()]);
	}
}

void DEventWriterROOT::Fill_ParticleData(bool locKinFitFlag, TTree* locTree, string locParticleBranchName, const DKinematicData* locKinematicData, const DKinematicData* locKinematicData_Measured, const DEventRFBunch* locEventRFBunch, const map<const DNeutralShower*, int>& locShowerToIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DDetectorMatches* locDetectorMatches) const
{
	//KINEMATICS: MEASURED
	TLorentzVector locX4_Measured(locKinematicData_Measured->position().X(), locKinematicData_Measured->position().Y(), locKinematicData_Measured->position().Z(), locKinematicData_Measured->time());
	Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4_Measured", &locX4_Measured, (*dTObjectMap)[locTree->GetName()]);
	TLorentzVector locP4_Measured(locKinematicData_Measured->momentum().X(), locKinematicData_Measured->momentum().Y(), locKinematicData_Measured->momentum().Z(), locKinematicData_Measured->energy());
	Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4_Measured", &locP4_Measured, (*dTObjectMap)[locTree->GetName()]);

	//KINEMATICS: KINFIT
	if(locKinFitFlag)
	{
		TLorentzVector locX4_KinFit(locKinematicData->position().X(), locKinematicData->position().Y(), locKinematicData->position().Z(), locKinematicData->time());
		Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "X4_KinFit", &locX4_KinFit, (*dTObjectMap)[locTree->GetName()]);
		TLorentzVector locP4_KinFit(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());
		Fill_TObjectData<TLorentzVector>(locTree, locParticleBranchName, "P4_KinFit", &locP4_KinFit, (*dTObjectMap)[locTree->GetName()]);
	}

	const DChargedTrackHypothesis* locChargedTrackHypothesis = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData);
	if(locChargedTrackHypothesis != NULL)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis_Measured = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData_Measured);

		//associated objects
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);

		const DBCALShower* locBCALShower = NULL;
		if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
			locBCALShower = dynamic_cast<const DBCALShower*>(locChargedTrackHypothesis->dBCALShowerMatchParams.dShowerObject);

		const DFCALShower* locFCALShower = NULL;
		if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
			locFCALShower = dynamic_cast<const DFCALShower*>(locChargedTrackHypothesis->dFCALShowerMatchParams.dShowerObject);

		//IDENTIFIER / MATCHING
		Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ObjectID", locChargedTrackHypothesis->candidateid);
		if(locMCThrownMatching != NULL)
		{
			Int_t locMatchID = -1;
			double locMatchFOM = 0.0;
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis_Measured, locMatchFOM);
			if((locMCThrown != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
				locMatchID = locThrownObjectIDMap.find(locMCThrown)->second;
			Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "MatchID", locMatchID);
		}

		//PID QUALITY
		Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Tracking", locChargedTrackHypothesis->dNDF_Track);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_Tracking", locChargedTrackHypothesis->dChiSq_Track);
		Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locChargedTrackHypothesis->dNDF_Timing);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locChargedTrackHypothesis_Measured->dChiSq_Timing);
		if(locKinFitFlag)
			Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locChargedTrackHypothesis->dChiSq_Timing);
		Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_DCdEdx", locChargedTrackHypothesis->dNDF_DCdEdx);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_DCdEdx", locChargedTrackHypothesis->dChiSq_DCdEdx);
		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "AvgBeta_Timing", locBeta_Timing);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "HitTime", locChargedTrackHypothesis->t1());

		//DEPOSITED ENERGY
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "dEdx_CDC", locTrackTimeBased->ddEdx_CDC);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "dEdx_FDC", locTrackTimeBased->ddEdx_FDC);
		double locTOFdEdx = (locChargedTrackHypothesis->dTOFHitMatchParams.dTrackTimeBased != NULL) ? locChargedTrackHypothesis->dTOFHitMatchParams.dEdx : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "dEdx_TOF", locTOFdEdx);
		double locSCdEdx = (locChargedTrackHypothesis->dSCHitMatchParams.dTrackTimeBased != NULL) ? locChargedTrackHypothesis->dSCHitMatchParams.dEdx : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "dEdx_ST", locSCdEdx);
		double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "Energy_BCAL", locBCALEnergy);
		double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "Energy_FCAL", locFCALEnergy);

		//Track DOCA to Shower - BCAL
		double locDOCAToShower_BCAL = 999.0;
		if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
			locDOCAToShower_BCAL = locChargedTrackHypothesis->dBCALShowerMatchParams.dDOCAToShower;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "TrackDOCA_BCAL", locDOCAToShower_BCAL);

		//Track DOCA to Shower - FCAL
		double locDOCAToShower_FCAL = 999.0;
		if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
			locDOCAToShower_FCAL = locChargedTrackHypothesis->dFCALShowerMatchParams.dDOCAToShower;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "TrackDOCA_FCAL", locDOCAToShower_FCAL);
	}
	else
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData);
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis_Measured = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData_Measured);

		//associated objects
		const DNeutralShower* locNeutralShower = NULL;
		locNeutralParticleHypothesis->GetSingleT(locNeutralShower);
		const DBCALShower* locBCALShower = NULL;
		locNeutralShower->GetSingleT(locBCALShower);
		const DFCALShower* locFCALShower = NULL;
		locNeutralShower->GetSingleT(locFCALShower);

		//IDENTIFIER / MATCHING
		Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "ObjectID", (locShowerToIDMap.find(locNeutralShower))->second);
		if(locMCThrownMatching != NULL)
		{
			Int_t locMatchID = -1;
			double locMatchFOM = 0.0;
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis_Measured, locMatchFOM);
			if((locMCThrown != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
				locMatchID = locThrownObjectIDMap.find(locMCThrown)->second;
			Fill_FundamentalData<Int_t>(locTree, locParticleBranchName, "MatchID", locMatchID);
		}

		//PID QUALITY
		Fill_FundamentalData<UInt_t>(locTree, locParticleBranchName, "NDF_Timing", locNeutralParticleHypothesis->dNDF);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_Measured", locNeutralParticleHypothesis_Measured->dChiSq);
		if(locKinFitFlag)
			Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "ChiSq_Timing_KinFit", locNeutralParticleHypothesis->dChiSq);
		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "AvgBeta_Timing", locBeta_Timing);
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "HitTime", locNeutralParticleHypothesis->t1());

		//DEPOSITED ENERGY
		double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "Energy_BCAL", locBCALEnergy);
		double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "Energy_FCAL", locFCALEnergy);

		//Track DOCA to Shower - BCAL
		double locDistanceToNearestTrack_BCAL = 999.0;
		if(locBCALShower != NULL)
		{
			if(!locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locDistanceToNearestTrack_BCAL))
				locDistanceToNearestTrack_BCAL = 999.0;
			if(locDistanceToNearestTrack_BCAL > 999.0)
				locDistanceToNearestTrack_BCAL = 999.0;
		}
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "TrackDOCA_BCAL", locDistanceToNearestTrack_BCAL);

		//Track DOCA to Shower - FCAL
		double locDistanceToNearestTrack_FCAL = 999.0;
		if(locFCALShower != NULL)
		{
			if(!locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistanceToNearestTrack_FCAL))
				locDistanceToNearestTrack_FCAL = 999.0;
			if(locDistanceToNearestTrack_FCAL > 999.0)
				locDistanceToNearestTrack_FCAL = 999.0;
		}
		Fill_FundamentalData<Double_t>(locTree, locParticleBranchName, "TrackDOCA_FCAL", locDistanceToNearestTrack_FCAL);
	}
}

void DEventWriterROOT::Fill_UnusedParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch, const map<const DNeutralShower*, int>& locShowerToIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DDetectorMatches* locDetectorMatches) const
{
	//KINEMATICS: MEASURED
	TLorentzVector locX4_Measured(locKinematicData->position().X(), locKinematicData->position().Y(), locKinematicData->position().Z(), locKinematicData->time());
	Fill_ClonesData<TLorentzVector>(locTree, "Unused", "X4_Measured", &locX4_Measured, locArrayIndex, (*dClonesArrayMap)[locTree->GetName()]);
	TLorentzVector locP4_Measured(locKinematicData->momentum().X(), locKinematicData->momentum().Y(), locKinematicData->momentum().Z(), locKinematicData->energy());
	Fill_ClonesData<TLorentzVector>(locTree, "Unused", "P4_Measured", &locP4_Measured, locArrayIndex, (*dClonesArrayMap)[locTree->GetName()]);

	const DChargedTrackHypothesis* locChargedTrackHypothesis = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData);
	if(locChargedTrackHypothesis != NULL)
	{
		//associated objects
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);

		const DBCALShower* locBCALShower = NULL;
		if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
			locBCALShower = dynamic_cast<const DBCALShower*>(locChargedTrackHypothesis->dBCALShowerMatchParams.dShowerObject);

		const DFCALShower* locFCALShower = NULL;
		if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
			locFCALShower = dynamic_cast<const DFCALShower*>(locChargedTrackHypothesis->dFCALShowerMatchParams.dShowerObject);

		//IDENTIFIERS / MATCHING
		Fill_FundamentalData<Int_t>(locTree, "Unused", "ObjectID", locChargedTrackHypothesis->candidateid, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		if(locMCThrownMatching != NULL)
		{
			Int_t locMatchID = -1;
			double locMatchFOM = 0.0;
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrackHypothesis, locMatchFOM);
			if((locMCThrown != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
				locMatchID = locThrownObjectIDMap.find(locMCThrown)->second;
			Fill_FundamentalData<Int_t>(locTree, "Unused", "MatchID", locMatchID, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		}
		Fill_FundamentalData<Int_t>(locTree, "Unused", "PID", PDGtype(locKinematicData->PID()), locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//PID QUALITY
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_Tracking", locChargedTrackHypothesis->dNDF_Track, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_Tracking", locChargedTrackHypothesis->dChiSq_Track, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_Timing", locChargedTrackHypothesis->dNDF_Timing, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_Timing", locChargedTrackHypothesis->dChiSq_Timing, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_DCdEdx", locChargedTrackHypothesis->dNDF_DCdEdx, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_DCdEdx", locChargedTrackHypothesis->dChiSq_DCdEdx, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locChargedTrackHypothesis, locEventRFBunch, true);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "AvgBeta_Timing", locBeta_Timing, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "HitTime", locChargedTrackHypothesis->t1(), locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//DEPOSITED ENERGY
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_CDC", locTrackTimeBased->ddEdx_CDC, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_FDC", locTrackTimeBased->ddEdx_FDC, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locTOFdEdx = (locChargedTrackHypothesis->dTOFHitMatchParams.dTrackTimeBased != NULL) ? locChargedTrackHypothesis->dTOFHitMatchParams.dEdx : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_TOF", locTOFdEdx, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locSCdEdx = (locChargedTrackHypothesis->dSCHitMatchParams.dTrackTimeBased != NULL) ? locChargedTrackHypothesis->dSCHitMatchParams.dEdx : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_ST", locSCdEdx, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "Energy_BCAL", locBCALEnergy, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "Energy_FCAL", locFCALEnergy, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//Track DOCA to Shower - BCAL
		double locDOCAToShower_BCAL = 999.0;
		if(locChargedTrackHypothesis->dBCALShowerMatchParams.dTrackTimeBased != NULL)
			locDOCAToShower_BCAL = locChargedTrackHypothesis->dBCALShowerMatchParams.dDOCAToShower;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "TrackDOCA_BCAL", locDOCAToShower_BCAL, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//Track DOCA to Shower - FCAL
		double locDOCAToShower_FCAL = 999.0;
		if(locChargedTrackHypothesis->dFCALShowerMatchParams.dTrackTimeBased != NULL)
			locDOCAToShower_FCAL = locChargedTrackHypothesis->dFCALShowerMatchParams.dDOCAToShower;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "TrackDOCA_FCAL", locDOCAToShower_FCAL, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
	}
	else
	{
		const DNeutralParticleHypothesis* locNeutralParticleHypothesis = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData);

		//associated objects
		const DNeutralShower* locNeutralShower = NULL;
		locNeutralParticleHypothesis->GetSingleT(locNeutralShower);
		const DBCALShower* locBCALShower = NULL;
		locNeutralShower->GetSingleT(locBCALShower);
		const DFCALShower* locFCALShower = NULL;
		locNeutralShower->GetSingleT(locFCALShower);

		//IDENTIFIERS
		Fill_FundamentalData<Int_t>(locTree, "Unused", "ObjectID", (locShowerToIDMap.find(locNeutralShower))->second, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		if(locMCThrownMatching != NULL)
		{
			Int_t locMatchID = -1;
			double locMatchFOM = 0.0;
			const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralParticleHypothesis, locMatchFOM);
			if((locMCThrown != NULL) && (locMatchFOM >= locMinThrownMatchFOM))
				locMatchID = locThrownObjectIDMap.find(locMCThrown)->second;
			Fill_FundamentalData<Int_t>(locTree, "Unused", "MatchID", locMatchID, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		}
		Fill_FundamentalData<Int_t>(locTree, "Unused", "PID", PDGtype(locKinematicData->PID()), locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//PID QUALITY
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_Tracking", 0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_Tracking", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_Timing", locNeutralParticleHypothesis->dNDF, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_Timing", locNeutralParticleHypothesis->dChiSq, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<UInt_t>(locTree, "Unused", "NDF_DCdEdx", 0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "ChiSq_DCdEdx", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locBeta_Timing = dAnalysisUtilities->Calc_Beta_Timing(locNeutralParticleHypothesis, locEventRFBunch);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "AvgBeta_Timing", locBeta_Timing, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "HitTime", locNeutralParticleHypothesis->t1(), locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//DEPOSITED ENERGY
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_CDC", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_FDC", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_TOF", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		Fill_FundamentalData<Double_t>(locTree, "Unused", "dEdx_ST", 0.0, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locBCALEnergy = (locBCALShower != NULL) ? locBCALShower->E : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "Energy_BCAL", locBCALEnergy, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
		double locFCALEnergy = (locFCALShower != NULL) ? locFCALShower->getEnergy() : 0.0;
		Fill_FundamentalData<Double_t>(locTree, "Unused", "Energy_FCAL", locFCALEnergy, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//Track DOCA to Shower - BCAL
		double locDistanceToNearestTrack_BCAL = 999.0;
		if(locBCALShower != NULL)
		{
			if(!locDetectorMatches->Get_DistanceToNearestTrack(locBCALShower, locDistanceToNearestTrack_BCAL))
				locDistanceToNearestTrack_BCAL = 999.0;
			if(locDistanceToNearestTrack_BCAL > 999.0)
				locDistanceToNearestTrack_BCAL = 999.0;
		}
		Fill_FundamentalData<Double_t>(locTree, "Unused", "TrackDOCA_BCAL", locDistanceToNearestTrack_BCAL, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);

		//Track DOCA to Shower - FCAL
		double locDistanceToNearestTrack_FCAL = 999.0;
		if(locFCALShower != NULL)
		{
			if(!locDetectorMatches->Get_DistanceToNearestTrack(locFCALShower, locDistanceToNearestTrack_FCAL))
				locDistanceToNearestTrack_FCAL = 999.0;
			if(locDistanceToNearestTrack_FCAL > 999.0)
				locDistanceToNearestTrack_FCAL = 999.0;
		}
		Fill_FundamentalData<Double_t>(locTree, "Unused", "TrackDOCA_FCAL", locDistanceToNearestTrack_FCAL, locArrayIndex, locMinArraySize, (*dNumUnusedArraySizeMap)[locTree]);
	}
}

