#include "DTreeInterface.h"


//Need one period: shared amongst threads
	//In plugin: class member
	//In action: static global? (map of unique-string to interface)
		//constructor of action would have to register threads. doable, but ugh.
	//In action: class member?
		//if branches already created, 
		//array sizes must be static global
		//num-threads-register must also be static global
		//static-global: 

int& DEventWriterROOT::Get_NumEventWriterThreads(void) const
{
	// must be read/used entirely in global root lock: when 0, close all files, modifying files needs global root lock
	static int locNumEventWriterThreads = 0;
	return locNumEventWriterThreads;
}

map<TTree*, map<string, unsigned int> >& DEventWriterROOT::Get_FundamentalArraySizeMap(void) const
{
	// all the trees (and thus these maps) are all created at once, within a global root lock
	// thus, don't need a lock when reading EITHER the inner or outer maps
	// however, when filling the tree the value will change, so modify within a file-lock
	static map<TTree*, map<string, unsigned int> > locFundamentalArraySizeMap;
	return locFundamentalArraySizeMap;
}

/********************************************************************* INITIALIZE *********************************************************************/

//Constructor
DTreeInterface::DTreeInterface(string locTreeName, string locFileName) : dFileName(locFileName)
{
	japp->RootWriteLock();
	{
		++Get_NumEventWriterThreads();
	}
	japp->RootUnLock();
	GetOrCreate_FileAndTree(locTreeName);
}

DTreeInterface::~DTreeInterface(void)
{
	japp->RootWriteLock();
	{
		--Get_NumEventWriterThreads();
		if(Get_NumEventWriterThreads() != 0)
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

void DEventWriterROOT::GetOrCreate_FileAndTree(string locTreeName) const
{
	japp->RootWriteLock();
	{
		//see if root file exists already
		TFile* locOutputFile = (TFile*)gROOT->GetListOfFiles()->FindObject(dFileName.c_str());
		if(locOutputFile == NULL)
			locOutputFile = new TFile(dFileName.c_str(), "RECREATE");
		else
			locOutputFile->cd();

		//see if ttree exists already
		TTree* locTree = (TTree*)gDirectory->Get(locTreeName.c_str());
		if(locTree != NULL)
		{
			dTree = locTree;
			japp->RootUnLock();
			return; //already created
		}

		//create tree
		dTree = new TTree(locTreeName.c_str(), locTreeName.c_str());
	}
	japp->RootUnLock();
}

void Create_Branches(const DTreeBranchRegister& locTreeBranchRegister)
{
	japp->WriteLock(dFileName); //LOCK FILE
	{
		//if there are branches already, don't do anything
		if(!dTree->GetListOfBranches()->IsEmpty())
		{
			japp->Unlock(dFileName); //UNLOCK FILE
			return;
		}

		//loop over branches
		for(auto locBranchIterator : locTreeFillData->dBranchTypeMap)
		{
			string locBranchName = locBranchIterator->first;
			type_index locTypeIndex = locBranchIterator->second;

			//array size
			auto locSizeIterator = locTreeFillData->dInitialArraySizeMap.find(locBranchName);
			if(locSizeIterator == locTreeFillData->dInitialArraySizeMap.end())
			{
				//not an array
				Create_Branch(locBranchName, locTypeIndex, 0, "");
				continue;
			}
			size_t locArraySize = locSizeIterator->second;

			//array size name
			auto locSizeNameIterator = locTreeFillData->dArraySizeNameMap.find(locBranchName);
			if(locSizeNameIterator == locTreeFillData->dArraySizeNameMap.end())
				Create_Branch(locBranchName, locTypeIndex, locArraySize, ""); //clones array
			else
				Create_Branch(locBranchName, locTypeIndex, locArraySize, locSizeNameIterator->second); //fundamental array
		}
	}
	japp->Unlock(dFileName); //UNLOCK FILE
}

/**************************************************************** FILL BRANCHES & TREE ****************************************************************/

void Fill(const DTreeFillData* locTreeFillData)
{
	japp->WriteLock(dFileName); //LOCK FILE
	{
		//loop over branches
		for(auto locBranchIterator : locTreeFillData->dFillData)
		{
			//type
			string locBranchName = locBranchIterator->first;
			type_index locTypeIndex = locBranchIterator->second.first;
			void* locVoidPointer = locBranchIterator->second.second;

			if(dTree->GetBranch(locBranchName.c_str()) == NULL)
			{
				cout << "WARNING, CANNOT FILL DATA, BRANCH " << locBranchName << " DOES NOT EXIST." << endl;
				continue;
			}

			//fill
			auto locNumFilledIterator = locTreeFillData->dArrayNumFilledMap.find(locBranchName);
			bool locIsArrayFlag = (locNumFilledIterator != locTreeFillData->dArrayNumFilledMap.end());
			if(locIsArrayFlag) //true: is array
			{
				//increase size if necessary
				unsigned int locCurrentArraySize = dFundamentalArraySizeMap[locBranchName]; //may not be in map! (tobj)
				size_t locNumFilled = locNumFilledIterator->second;
				if(locNumFilled >= locCurrentArraySize)
					Increase_ArraySize(locBranchName, locTypeIndex, locNumFilled)

				for(size_t locArrayIndex = 0; locArrayIndex <= locNumFilledIterator->second; ++locArrayIndex)
					Fill(locBranchName, locTypeIndex, locVoidPointer, locIsArrayFlag, locArrayIndex);
			}
			else // is not an array
				Fill(locBranchName, locTypeIndex, locVoidPointer, false);
		}

		//fill tree
		dTree->Fill();
	}
	japp->Unlock(dFileName); //UNLOCK FILE
}

void Increase_ArraySize(string locBranchName, type_index locTypeIndex, size_t locNewArraySize)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(const char*)))
		Increase_ArraySize<const char*>(locBranchName, locNewArraySize);
	else if(locTypeIndex == type_index(typeid(Char_t)))
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

void Fill(string locBranchName, type_index locTypeIndex, void* locVoidPointer, bool locIsArrayFlag, size_t locArrayIndex)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(const char*)))
		Get_Pointer_Fundamental<const char*>(locBranchName)[locArrayIndex] = *(static_cast<const char*>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Char_t)))
		Get_Pointer_Fundamental<Char_t>(locBranchName)[locArrayIndex] = *(static_cast<Char_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Get_Pointer_Fundamental<UChar_t>(locBranchName)[locArrayIndex] = *(static_cast<UChar_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Get_Pointer_Fundamental<Short_t>(locBranchName)[locArrayIndex] = *(static_cast<Short_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Get_Pointer_Fundamental<UShort_t>(locBranchName)[locArrayIndex] = *(static_cast<UShort_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Get_Pointer_Fundamental<Int_t>(locBranchName)[locArrayIndex] = *(static_cast<Int_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Get_Pointer_Fundamental<UInt_t>(locBranchName)[locArrayIndex] = *(static_cast<UInt_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Get_Pointer_Fundamental<Float_t>(locBranchName)[locArrayIndex] = *(static_cast<Float_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Get_Pointer_Fundamental<Double_t>(locBranchName)[locArrayIndex] = *(static_cast<Double_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Get_Pointer_Fundamental<Long64_t>(locBranchName)[locArrayIndex] = *(static_cast<Long64_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Get_Pointer_Fundamental<ULong64_t>(locBranchName)[locArrayIndex] = *(static_cast<ULong64_t>(locVoidPointer));
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Get_Pointer_Fundamental<Bool_t>(locBranchName)[locArrayIndex] = *(static_cast<Bool_t>(locVoidPointer));

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

void Create_Branch(string locBranchName, type_index locTypeIndex, size_t locArraySize, string locArraySizeName)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(const char*)))
		Create_Branch<const char*>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Char_t)))
		Create_Branch<Char_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Create_Branch<UChar_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Create_Branch<Short_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Create_Branch<UShort_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Create_Branch<Int_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Create_Branch<UInt_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Create_Branch<Float_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Create_Branch<Double_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Create_Branch<Long64_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Create_Branch<ULong64_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Create_Branch<Bool_t>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);

	//TObject
	else if(locTypeIndex == type_index(typeid(TVector2)))
		Create_Branch<TVector2>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(TVector3)))
		Create_Branch<TVector3>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
	else if(locTypeIndex == type_index(typeid(TLorentzVector)))
		Create_Branch<TLorentzVector>(locBranchName, locTypeIndex, locArraySize, locArraySizeName);
}




		template <typename DType> void Create_Branch_Fundamental(string locBranchName);
		template <typename DType> void Create_Branch_TObject(string locBranchName);
		template <typename DType> void Create_Branch_FundamentalArray(string locBranchName, string locArraySizeString, unsigned int locInitialSize);
		template <typename DType> void Create_Branch_ClonesArray(string locBranchName, unsigned int locSize);


