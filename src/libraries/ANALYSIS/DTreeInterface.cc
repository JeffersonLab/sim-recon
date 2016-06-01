#include "DTreeInterface.h"

/************************************************* STATIC-VARIABLE-ACCESSING PRIVATE MEMBER FUNCTIONS *************************************************/

map<string, int>& DTreeInterface::Get_NumWritersByFileMap(void) const
{
	// string is file name
	// must be read/used entirely in global root lock: when 0, close the file (modifying files needs global root lock)
	static map<string, int> locNumWritersByFileMap;
	return locNumWritersByFileMap;
}

map<string, size_t>& DTreeInterface::Get_FundamentalArraySizeMap(TTree* locTree) const
{
	//A "global" lock (shared across all threads & DTreeInterface objects) is necessary to read the outer map
		//This lock is acquired locally: DO NOT CALL THIS FUNCTION WHILE WITHIN A LOCK!!
	//However, only a file-lock is necessary to read the inner map once it is acquired:
		//Only the given tree is viewing/modifying it

	japp->WriteLock("DTreeInterface_SizeMap"); //LOCK MAP

	static map<TTree*, map<string, size_t> > locFundamentalArraySizeMap;
	map<string, size_t>& locTreeSpecificMap = locFundamentalArraySizeMap[locTree];

	japp->Unlock("DTreeInterface_SizeMap"); //UNLOCK MAP

	return locTreeSpecificMap;
}

/********************************************************************* INITIALIZE *********************************************************************/

DTreeInterface* DTreeInterface::Create_DTreeInterface(string locTreeName, string locFileName)
{
	if(locFileName == "hd_root.root")
	{
		cout << "WARNING: SAVING TREES TO hd_root.root IS NOT SUPPORTED (or wise). ";
		cout << "RETURNING NULL FROM DTreeInterface::Create_DTreeInterface()" << endl;
		return NULL;
	}
	return new DTreeInterface(locTreeName, locFileName);
}

//Constructor
DTreeInterface::DTreeInterface(string locTreeName, string locFileName) : dFileName(locFileName)
{
	japp->RootWriteLock();
	{
		map<string, int>& locNumWritersByFileMap = Get_NumWritersByFileMap();
		if(locNumWritersByFileMap.find(dFileName) == locNumWritersByFileMap.end())
			locNumWritersByFileMap[dFileName] = 1;
		else
			++locNumWritersByFileMap[dFileName];
	}
	japp->RootUnLock();

	GetOrCreate_FileAndTree(locTreeName);
}

//Destructor
DTreeInterface::~DTreeInterface(void)
{
	japp->RootWriteLock();
	{
		map<string, int>& locNumWritersByFileMap = Get_NumWritersByFileMap();
		--locNumWritersByFileMap[dFileName];
		if(locNumWritersByFileMap[dFileName] != 0)
		{
			japp->RootUnLock();
			return;
		}

		//see if root file exists already
		TFile* locOutputFile = (TFile*)gROOT->GetListOfFiles()->FindObject(dFileName.c_str());
		locOutputFile->Write(0, TObject::kOverwrite);
		locOutputFile->Close();
		delete locOutputFile;
	}
	japp->RootUnLock();
}

void DTreeInterface::GetOrCreate_FileAndTree(string locTreeName)
{
	japp->RootWriteLock();
	{
		TDirectory* locCurrentDir = gDirectory;

		//see if root file exists already
		TFile* locOutputFile = (TFile*)gROOT->GetListOfFiles()->FindObject(dFileName.c_str());
		if(locOutputFile == nullptr)
			locOutputFile = new TFile(dFileName.c_str(), "RECREATE");
		locOutputFile->cd();

		//see if ttree exists already. if not, create it
		dTree = (TTree*)gDirectory->Get(locTreeName.c_str());
		if(dTree == nullptr)
			dTree = new TTree(locTreeName.c_str(), locTreeName.c_str());

		locCurrentDir->cd();
	}
	japp->RootUnLock();
}

/****************************************************************** CREATE BRANCHES *******************************************************************/

bool DTreeInterface::Get_BranchesCreatedFlag(void) const
{
	japp->WriteLock(dFileName); //LOCK FILE
	bool locBranchesCreatedFlag = (dTree->GetNbranches() > 0);
	japp->Unlock(dFileName); //UNLOCK FILE
	return locBranchesCreatedFlag;
}

const TList* DTreeInterface::Get_UserInfo(void) const
{
	if(Get_BranchesCreatedFlag())
		return dTree->GetUserInfo();

	cout << "WARNING: CANNOT GET USER INFO BEFORE BRANCHES CREATED. RETURNING NULL IN DTreeInterface::Get_UserInfo()" << endl;
	return NULL; //NOT SUPPORTED! //Unsafe otherwise. This guarantees that the user info is setup first, and won't be modified while reading it
}

void DTreeInterface::Create_Branches(const DTreeBranchRegister& locTreeBranchRegister)
{
	//MUST CARRY AROUND A REFERENCE TO THIS.  ONLY READ/MODIFY THE MAP WITHIN A FILE LOCK. 
	map<string, size_t>& locFundamentalArraySizeMap = Get_FundamentalArraySizeMap(dTree);

	japp->WriteLock(dFileName); //LOCK FILE
	{
		//if there are branches already, don't do anything
		if(dTree->GetNbranches() > 0)
		{
			japp->Unlock(dFileName); //UNLOCK FILE
			return;
		}

		//set user info
		TList* locInputUserInfo = locTreeBranchRegister.Get_UserInfo();
		TList* locTreeUserInfo = dTree->GetUserInfo();
		for(Int_t loc_i = 0; loc_i < locInputUserInfo->GetSize(); ++loc_i)
			locTreeUserInfo->Add(locInputUserInfo->At(loc_i));

		//loop over branches
			//branches that are arrays, for which the branch with the size was not created it yet
			//save them, and create them at the end
		vector<string> locPostponedBranchNames;
		for(const auto& locBranchName : locTreeBranchRegister.dBranchNames)
		{
			//check if need to postpone this branch
			auto locSizeNameIterator = locTreeBranchRegister.dArraySizeNameMap.find(locBranchName);
			if(locSizeNameIterator != locTreeBranchRegister.dArraySizeNameMap.end())
			{
				//fundamental array: make sure array size branch has been created first. if not, postpone it
				if(dTree->GetBranch(locSizeNameIterator->second.c_str()) == nullptr)
				{
					locPostponedBranchNames.push_back(locBranchName);
					continue;
				}
			}

			//it's ok: create it
			Create_Branch(locTreeBranchRegister, locBranchName, locFundamentalArraySizeMap);
		}

		//create postponed branches
		for(string& locBranchName : locPostponedBranchNames)
			Create_Branch(locTreeBranchRegister, locBranchName, locFundamentalArraySizeMap);
	}
	japp->Unlock(dFileName); //UNLOCK FILE
}

void DTreeInterface::Create_Branch(const DTreeBranchRegister& locTreeBranchRegister, string locBranchName, map<string, size_t>& locFundamentalArraySizeMap)
{
	const type_index& locTypeIndex = locTreeBranchRegister.dBranchTypeMap.find(locBranchName)->second;

	//array size
	auto locSizeIterator = locTreeBranchRegister.dInitialArraySizeMap.find(locBranchName);
	if(locSizeIterator == locTreeBranchRegister.dInitialArraySizeMap.end())
	{
		//not an array
		Create_Branch(locBranchName, locTypeIndex, 0, "");
		return;
	}
	size_t locArraySize = locSizeIterator->second;

	//array size name
	auto locSizeNameIterator = locTreeBranchRegister.dArraySizeNameMap.find(locBranchName);
	if(locSizeNameIterator == locTreeBranchRegister.dArraySizeNameMap.end())
		Create_Branch(locBranchName, locTypeIndex, locArraySize, ""); //clones array
	else //fundamental array
	{
		Create_Branch(locBranchName, locTypeIndex, locArraySize, locSizeNameIterator->second);
		locFundamentalArraySizeMap[locBranchName] = locArraySize;
	}
}

void DTreeInterface::Create_Branch(string locBranchName, type_index locTypeIndex, size_t locArraySize, string locArraySizeName)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(Char_t)))
		Create_Branch<Char_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Create_Branch<UChar_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Create_Branch<Short_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Create_Branch<UShort_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Create_Branch<Int_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Create_Branch<UInt_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Create_Branch<Float_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Create_Branch<Double_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Create_Branch<Long64_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Create_Branch<ULong64_t>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Create_Branch<Bool_t>(locBranchName, locArraySize, locArraySizeName);

	//TObject
	else if(locTypeIndex == type_index(typeid(TVector2)))
		Create_Branch<TVector2>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(TVector3)))
		Create_Branch<TVector3>(locBranchName, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(TLorentzVector)))
		Create_Branch<TLorentzVector>(locBranchName, locArraySize, locArraySizeName);
}

/**************************************************************** FILL BRANCHES & TREE ****************************************************************/

void DTreeInterface::Fill(DTreeFillData& locTreeFillData)
{
	//MUST CARRY AROUND A REFERENCE TO THIS.  ONLY READ/MODIFY THE MAP WITHIN A FILE LOCK. 
	map<string, size_t>& locFundamentalArraySizeMap = Get_FundamentalArraySizeMap(dTree);

	japp->WriteLock(dFileName); //LOCK FILE
	{
		//loop over branches
		for(auto locBranchPair : *(locTreeFillData.dFillData))
		{
			//type
			string locBranchName = locBranchPair.first;
			type_index locTypeIndex = locBranchPair.second.first;
			deque<void*>& locVoidDeque = *(locBranchPair.second.second);

			if(dTree->GetBranch(locBranchName.c_str()) == NULL)
			{
				cout << "WARNING, CANNOT FILL DATA, BRANCH " << locBranchName << " DOES NOT EXIST." << endl;
				continue;
			}

			//check if is array. if not, fill
			auto locLargestIndexFilledIterator = locTreeFillData.dArrayLargestIndexFilledMap->find(locBranchName);
			bool locIsArrayFlag = (locLargestIndexFilledIterator != locTreeFillData.dArrayLargestIndexFilledMap->end());
			if(!locIsArrayFlag)
			{
				Fill(locBranchName, locTypeIndex, locVoidDeque[0], false);
				continue;
			}

			//is array, get how many to fill
			size_t& locLargestIndexFilled = locLargestIndexFilledIterator->second;

			//increase array size if necessary
			auto locFundamentalArraySizeIterator = locFundamentalArraySizeMap.find(locBranchName);
			if(locFundamentalArraySizeIterator != locFundamentalArraySizeMap.end())
			{
				size_t locCurrentArraySize = locFundamentalArraySizeMap[locBranchName]; //may not be in map! (tobj)
				if((locLargestIndexFilled + 1) >= locCurrentArraySize)
					Increase_ArraySize(locBranchName, locTypeIndex, locLargestIndexFilled + 1);
				locFundamentalArraySizeMap[locBranchName] = locLargestIndexFilled + 1;
			}
			else //is clones array: clear it
				Get_Pointer_TClonesArray(locBranchName)->Clear(); //empties array

			//fill array
			for(size_t locArrayIndex = 0; locArrayIndex <= locLargestIndexFilled; ++locArrayIndex)
				Fill(locBranchName, locTypeIndex, locVoidDeque[locArrayIndex], true, locArrayIndex);

			//reset DTreeFillData for next event!
			locLargestIndexFilled = 0;
cout << "value = " << (*locTreeFillData.dArrayLargestIndexFilledMap)[locBranchName] << endl;
		}

		//fill tree
		dTree->Fill();
	}
	japp->Unlock(dFileName); //UNLOCK FILE
}

void DTreeInterface::Increase_ArraySize(string locBranchName, type_index locTypeIndex, size_t locNewArraySize)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(Char_t)))
		Increase_ArraySize<Char_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Increase_ArraySize<UChar_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Increase_ArraySize<Short_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Increase_ArraySize<UShort_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Increase_ArraySize<Int_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Increase_ArraySize<UInt_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Increase_ArraySize<Float_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Increase_ArraySize<Double_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Increase_ArraySize<Long64_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Increase_ArraySize<ULong64_t>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Increase_ArraySize<Bool_t>(locBranchName, locNewArraySize);
}

void DTreeInterface::Fill(string locBranchName, type_index locTypeIndex, void* locVoidPointer, bool locIsArrayFlag, size_t locArrayIndex)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(Char_t)))
		Get_Pointer_Fundamental<Char_t>(locBranchName)[locArrayIndex] = *(static_cast<Char_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Get_Pointer_Fundamental<UChar_t>(locBranchName)[locArrayIndex] = *(static_cast<UChar_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Get_Pointer_Fundamental<Short_t>(locBranchName)[locArrayIndex] = *(static_cast<Short_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Get_Pointer_Fundamental<UShort_t>(locBranchName)[locArrayIndex] = *(static_cast<UShort_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Get_Pointer_Fundamental<Int_t>(locBranchName)[locArrayIndex] = *(static_cast<Int_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Get_Pointer_Fundamental<UInt_t>(locBranchName)[locArrayIndex] = *(static_cast<UInt_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Get_Pointer_Fundamental<Float_t>(locBranchName)[locArrayIndex] = *(static_cast<Float_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Get_Pointer_Fundamental<Double_t>(locBranchName)[locArrayIndex] = *(static_cast<Double_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Get_Pointer_Fundamental<Long64_t>(locBranchName)[locArrayIndex] = *(static_cast<Long64_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Get_Pointer_Fundamental<ULong64_t>(locBranchName)[locArrayIndex] = *(static_cast<ULong64_t*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Get_Pointer_Fundamental<Bool_t>(locBranchName)[locArrayIndex] = *(static_cast<Bool_t*>(locVoidPointer));

	//TObject
	else if(locTypeIndex == type_index(typeid(TVector3)))
	{
		TVector3& locObject = *(static_cast<TVector3*>(locVoidPointer));
		Fill_TObject<TVector3>(locBranchName, locObject, locIsArrayFlag, locArrayIndex);
	}
	else if(locTypeIndex == type_index(typeid(TVector2)))
	{
		TVector2& locObject = *(static_cast<TVector2*>(locVoidPointer));
		Fill_TObject<TVector2>(locBranchName, locObject, locIsArrayFlag, locArrayIndex);
	}
	else if(locTypeIndex == type_index(typeid(TLorentzVector)))
	{
		TLorentzVector& locObject = *(static_cast<TLorentzVector*>(locVoidPointer));
		Fill_TObject<TLorentzVector>(locBranchName, locObject, locIsArrayFlag, locArrayIndex);
	}
}
