#ifndef DTreeInterface_h
#define DTreeInterface_h

#include <map>
#include <set>
#include <string>
#include <iostream>
#include <sstream>

#include <TROOT.h>
#include <TTree.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <particleType.h>

using namespace std;


//Need one period: shared amongst threads
	//In plugin: class member
	//In action: static global? (map of unique-string to interface)
		//constructor of action would have to register threads. doable, but ugh.
	//In action: class member?
		//if branches already created, 
		//array sizes must be static global
		//num-threads-register must also be static global
		//static-global: 
class DTreeInterface
{
	//ASSUME: One leaf per branch: No splitting
	//ASSUME: Tree only stores: Fundamental type objects (not const char*!!), Fundamental type arrays, TObject's, TClonesArray's
	//ASSUME: TObject's are of type TVector3 & TLorentzVector: Need to expand template type calls if more are used

	public:

		/**************************************************************** INITIALIZE ****************************************************************/

		//Constructor
		DTreeInterface(string locTreeName, string locFileName);

		// Is optional. If needed array size is not specified, it will be set to a default value.
		void Set_InitialArraySize(string locArraySizeBranchName, UInt_t locInitialSize);

		void Create_Branches(const DTreeBranchRegister& locTreeBranchRegister);

		/************************************************************* GET ENTRY AND FILL ***********************************************************/

		void Fill(DTreeFillData* locTreeFillData);

		/************************************************************ GET BRANCHES AND DATA *********************************************************/

		//GET OBJECT POINTERS
		template <typename DType> DType* Get_Pointer_Fundamental(string locBranchName) const;
		template <typename DType> DType* Get_Pointer_TObject(string locBranchName) const;
		TClonesArray* Get_Pointer_TClonesArray(string locBranchName);

		/******************************************************** CREATE & FILL NEW BRANCHES ********************************************************/

		//CREATE BRANCHES
		template <typename DType> void Create_Branch_Fundamental(string locBranchName);
		template <typename DType> void Create_Branch_TObject(string locBranchName);
		template <typename DType> void Create_Branch_FundamentalArray(string locBranchName, string locArraySizeString, unsigned int locInitialSize);
		template <typename DType> void Create_Branch_ClonesArray(string locBranchName, unsigned int locSize);

		//FILL BRANCHES
		template <typename DType> void Fill_Fundamental(string locBranchName, DType locValue);
		template <typename DType> void Fill_Fundamental(string locBranchName, DType locValue, unsigned int locArrayIndex);
		template <typename DType> void Fill_TObject(string locBranchName, DType& locObject, unsigned int locArrayIndex);
		template <typename DType> void Fill_TObject(string locBranchName, DType& locObject);

	private:

		/************************************************************** MISCELLANEOUS ***************************************************************/

		DTreeInterface(void); //private default constructor: cannot call

		//For ROOT type string for fundamental data variables
			//Defined in https://root.cern.ch/root/htmldoc/TTree.html
		template<typename DType> struct DROOTTypeString { static const char* GetTypeString() {return "";} }; // Main template class

		/************************************************************* MEMBER VARIABLES *************************************************************/

		TTree* dTree;
		string dFileName;

		/****************************************************** BRANCH MEMORY AND TYPE MAPPING ******************************************************/

		//These objects are kept here because they must be kept somewhere: 
			//branches addresses of pointers, so pointers must reside somewhere permanent
			//However, for fundamental objects/arrays: memory stored in the branches themselves: don't need to hold onto them
		map<string, TClonesArray*> dMemoryMap_ClonesArray;
		map<string, TObject*> dMemoryMap_TObject;

		/************************************************************ ARRAY SIZE MAPPING ************************************************************/

		//Branch Fundamental Array Maps
		map<string, UInt_t> dFundamentalArraySizeMap; //keys are the names of the branches that contain sizes, value is the current size //make sure has init value!!
};

/**************************************************************** GET BRANCHES AND DATA ***************************************************************/

//GET POINTERS
template <typename DType> inline DType* DTreeInterface::Get_Pointer_Fundamental(string locBranchName) const
{
	TBranch* locBranch = dTree->GetBranch(locBranchName.c_str());
	return ((locBranch != NULL) ? (DType*)locBranch->GetAddress() : NULL);
}

template <typename DType> inline DType* DTreeInterface::Get_Pointer_TObject(string locBranchName) const
{
	TBranch* locBranch = dTree->GetBranch(locBranchName.c_str());
	return ((locBranch != NULL) ? *(DType**)locBranch->GetAddress() : NULL);
}

inline TClonesArray* DTreeInterface::Get_Pointer_TClonesArray(string locBranchName)
{
	TBranch* locBranch = dTree->GetBranch(locBranchName.c_str());
	return ((locBranch != NULL) ? *(TClonesArray**)locBranch->GetAddress() : NULL);
}

/************************************************************* CREATE & FILL NEW BRANCHES *************************************************************/

//CREATE BRANCHES
template <typename DType> inline void DTreeInterface::Create_Branch_Fundamental(string locBranchName)
{
	if(dTree->GetBranch(locBranchName.c_str()) != NULL)
		return; //already created

	string locTypeString = DROOTTypeString<DType>::GetTypeString();
	string locTypeName = locBranchName + string("/") + locTypeString;
	dTree->Branch(locBranchName.c_str(), new DType(), locTypeName.c_str());
}

template <typename DType> inline void DTreeInterface::Create_Branch_TObject(string locBranchName)
{
	if(dTree->GetBranch(locBranchName.c_str()) != NULL)
		return; //already created

	dMemoryMap_TObject[locBranchName] = (TObject*)(new DType());
	dTree->Branch(locBranchName.c_str(), (DType**)&(dMemoryMap_TObject[locBranchName]), 32000, 0); //0: don't split
}

template <typename DType> inline void DTreeInterface::Create_Branch_FundamentalArray(string locBranchName, string locArraySizeString, unsigned int locInitialSize)
{
	if(dTree->GetBranch(locBranchName.c_str()) != NULL)
		return; //already created

	string locTypeString = DROOTTypeString<DType>::GetTypeString();
	string locArrayName = locBranchName + string("[") + locArraySizeString + string("]/") + locTypeString;
	dTree->Branch(locBranchName.c_str(), new DType[locInitialSize], locArrayName.c_str());
	dFundamentalArraySizeMap[locBranchName] = locInitialSize;
}

template <typename DType> inline void DTreeInterface::Create_Branch_ClonesArray(string locBranchName, unsigned int locSize)
{
	if(dTree->GetBranch(locBranchName.c_str()) != NULL)
		return; //already created

	dMemoryMap_ClonesArray[locBranchName] = new TClonesArray(DType::Class()->GetName(), locSize);
	dTree->Branch(locBranchName.c_str(), &(dMemoryMap_ClonesArray[locBranchName]), 32000, 0); //0: don't split
}

//FILL BRANCHES
template <typename DType> inline void DTreeInterface::Fill_Fundamental(string locBranchName, DType locValue)
{
	*Get_Pointer_Fundamental<DType>(locBranchName) = locValue;
}

template <typename DType> inline void DTreeInterface::Fill_Fundamental(string locBranchName, DType locValue, unsigned int locArrayIndex)
{
	//create a new, larger array if the current one is too small
		//would rather do this in advance, but don't assume that the user has done so!
	unsigned int locCurrentArraySize = dFundamentalArraySizeMap[locBranchName];
	DType* locArray = Get_Pointer_Fundamental<DType>(locBranchName);
	if(locArrayIndex >= locCurrentArraySize)
	{
		DType* locOldArray = locArray;
		dTree->SetBranchAddress(locBranchName.c_str(), new DType[locArrayIndex + 1]);
		locArray = Get_Pointer_Fundamental<DType>(locBranchName);

		//copy the old contents into the new array
		for(unsigned int loc_i = 0; loc_i < locCurrentArraySize; ++loc_i)
			locArray[loc_i] = locOldArray[loc_i];

		delete[] locOldArray;
		dFundamentalArraySizeMap[locBranchName] = locArrayIndex + 1;
	}

	//set the data
	locArray[locArrayIndex] = locValue;
}

template <typename DType> inline void DTreeInterface::Fill_TObject(string locBranchName, DType& locObject, unsigned int locArrayIndex)
{
	TClonesArray* locClonesArray = Get_Pointer_TClonesArray(locBranchName);
	if(locArrayIndex == 0)
		locClonesArray->Clear(); //empties array
	*(DType*)locClonesArray->ConstructedAt(locArrayIndex) = locObject;
}

template <typename DType> inline void DTreeInterface::Fill_TObject(string locBranchName, DType& locObject)
{
	*Get_Pointer_TObject<DType>(locBranchName) = locObject;
}

/************************************************************** TEMPLATE SPECIALIZATIONS **************************************************************/

template<> struct DTreeInterface::DROOTTypeString<const char*> { static const char* GetTypeString() {return "C";} };
template<> struct DTreeInterface::DROOTTypeString<Char_t> { static const char* GetTypeString() {return "B";} };
template<> struct DTreeInterface::DROOTTypeString<UChar_t> { static const char* GetTypeString() {return "b";} };
template<> struct DTreeInterface::DROOTTypeString<Short_t> { static const char* GetTypeString() {return "S";} };
template<> struct DTreeInterface::DROOTTypeString<UShort_t> { static const char* GetTypeString() {return "s";} };
template<> struct DTreeInterface::DROOTTypeString<Int_t> { static const char* GetTypeString() {return "I";} };
template<> struct DTreeInterface::DROOTTypeString<UInt_t> { static const char* GetTypeString() {return "i";} };
template<> struct DTreeInterface::DROOTTypeString<Float_t> { static const char* GetTypeString() {return "F";} };
template<> struct DTreeInterface::DROOTTypeString<Double_t> { static const char* GetTypeString() {return "D";} };
template<> struct DTreeInterface::DROOTTypeString<Long64_t> { static const char* GetTypeString() {return "L";} };
template<> struct DTreeInterface::DROOTTypeString<ULong64_t> { static const char* GetTypeString() {return "l";} };
template<> struct DTreeInterface::DROOTTypeString<Bool_t> { static const char* GetTypeString() {return "O";} };


//INCREASE ARRAY SIZE //For when reading only! Not when writing
template <typename DType> inline void DTreeInterface::Increase_ArraySize(string locBranchName, int locNewArraySize)
{
	//create a new, larger array if the current one is too small
		//DOES NOT copy the old results!  In other words, only call BETWEEN entries, not DURING an entry
	DType* locOldBranchAddress = Get_Pointer_Fundamental<DType>(locBranchName);
	dTree->SetBranchAddress(locBranchName.c_str(), new DType[locNewArraySize]);
	delete[] locOldBranchAddress;
	dFundamentalArraySizeMap[locBranchName] = locNewArraySize;
}

template <typename DType> inline void Fill_TObject(string locBranchName, DType& locObject, bool locIsArrayFlag, size_t locArrayIndex)
{
	if(locIsArrayFlag)
		Fill_TObject<DType>(locBranchName, locObject, locArrayIndex);
	else
		Fill_TObject<DType>(locBranchName, locObject);
}

#endif //DTreeInterface_h

