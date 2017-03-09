#ifndef DTreeInterface_h
#define DTreeInterface_h

#include <map>
#include <typeindex>
#include <typeinfo>
#include <type_traits>
#include <string>
#include <iostream>
#include <sstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <JANA/JApplication.h>
#include <particleType.h>

#include "DTreeInterfaceObjects.h"

using namespace std;

class DTreeInterface
{
	//WARNING: This class ASSUMES that the only things you are saving to the given output files are trees managed by DTreeInterface objects.  
		//Saving anything else to these files will probably not work!  And would be unwise anyway ...

	//Usage: Create anywhere, as long as it will live for the life of the time you want to fill the tree. 
		//Regardless of how many threads act in that part of the code. 
		//Just delete it whenever you're done. TTree is saved and TFile is closed when all interfaces to a given file are deleted. 
		//E.g. Create in plugin/factory init(), delete in plugin/factory fini()
		//E.g. Create in analysis action Initialize(), delete in analysis action destructor

	//ASSUME: One leaf per branch: No splitting
	//ASSUME: Tree only stores: Fundamental type objects (not const char*!!), Fundamental type arrays, TObject's, TClonesArray's
	//ASSUME: TObject's are of type TVector3 & TLorentzVector: Need to expand template type calls if more are used

	public:

		/**************************************************************** INITIALIZE ****************************************************************/

		//Only public way to construct a DTreeInterface //Forces allocation on the heap
		//MUST DELETE WHEN FINISHED: OR ELSE DATA WON'T BE SAVED!!!
		static DTreeInterface* Create_DTreeInterface(string locTreeName, string locFileName);

		//Destructor
		~DTreeInterface(void);

		// Is optional. If needed array size is not specified, it will be set to a default value.
		void Set_InitialArraySize(string locArraySizeBranchName, UInt_t locInitialSize);

		bool Create_Branches(const DTreeBranchRegister& locTreeBranchRegister);

		void Set_TreeIndexBranchNames(string locTreeIndex_MajorBranchName, string locTreeIndex_MinorBranchName = "0");

		//Check/read info
		bool Get_BranchesCreatedFlag(void) const;
		const TList* Get_UserInfo(void) const;

		/******************************************************************* FILL *******************************************************************/

		void Fill(DTreeFillData& locTreeFillData); //not const: needs to reset arrays

	private:

		/**************************************************************** INITIALIZE ****************************************************************/

		//Constructors
		DTreeInterface(string locTreeName, string locFileName);
		DTreeInterface(void); //private default constructor: cannot call

		//Init
		void GetOrCreate_FileAndTree(string locTreeName);

		/*************************************************************** MISCELLANEOUS **************************************************************/

		//For ROOT type string for fundamental data variables
			//Defined in https://root.cern.ch/root/htmldoc/TTree.html
		template<typename DType> struct DROOTTypeString { static const char* GetTypeString() {return "";} }; // Main template class

		/************************************************************** CREATE BRANCHES *************************************************************/

		void Create_Branch(const DTreeBranchRegister& locTreeBranchRegister, string locBranchName, map<string, size_t>& locFundamentalArraySizeMap);
		void Create_Branch(string locBranchName, type_index locTypeIndex, size_t locArraySize, string locArraySizeName);

		//Enable this version if type inherits from TObject //void: is return type
		template <typename DType> typename enable_if<std::is_base_of<TObject, DType>::value, void>::type
				Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName);

		//Enable this version if type does NOT inherit from TObject //void: is return type
		template <typename DType> typename enable_if<!std::is_base_of<TObject, DType>::value, void>::type
				Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName);

		template <typename DType> void Create_Branch_Fundamental(string locBranchName);
		template <typename DType> void Create_Branch_TObject(string locBranchName);
		template <typename DType> void Create_Branch_FundamentalArray(string locBranchName, string locArraySizeString, unsigned int locInitialSize);
		template <typename DType> void Create_Branch_ClonesArray(string locBranchName, unsigned int locSize);

		/******************************************************************* FILL *******************************************************************/

		void Increase_ArraySize(string locBranchName, type_index locTypeIndex, size_t locNewArraySize);
		template <typename DType> void Increase_ArraySize(string locBranchName, int locNewArraySize);
		void Fill(string locBranchName, type_index locTypeIndex, void* locVoidPointer, bool locIsArrayFlag, size_t locArrayIndex = 0);
		template <typename DType> void Fill_TObject(string locBranchName, DType& locObject, bool locIsArrayFlag, size_t locArrayIndex);

		/*************************************************************** GET POINTERS ***************************************************************/

		template <typename DType> DType* Get_Pointer_Fundamental(string locBranchName) const;
		template <typename DType> DType* Get_Pointer_TObject(string locBranchName) const;
		TClonesArray* Get_Pointer_TClonesArray(string locBranchName);

		/******************************************** STATIC-VARIABLE-ACCESSING PRIVATE MEMBER FUNCTIONS ********************************************/

		//Some variables needs to be shared amongst threads (e.g. the array sizes for the branches)
		//However, you cannot make them global/extern/static/static-member variables in the header file:
			//They would be in the header file, and the header file is included in the ANALYSIS library AND in each plugin that uses it
				//When a header file is included in a src file, it's contents are essentially copied directly into it
			//Thus there are two instances of each static variable: one in each translation unit (library)
			//Supposedly(?) they are linked together during runtime when loading, so there is (supposedly) no undefined behavior.
			//However, this causes a double free (double-deletion) when these libraries are closed at the end of the program, crashing it.
		//Thus the variables must be in a single source file that is compiled into a single library
		//However, you (somehow?) cannot make them global/extern variables in the source file
			//This also (somehow?) causes the double-free problem above for (at least) stl containers
			//It works for pointers-to-stl-containers and fundamental types, but I dunno why.
			//It's not good encapsulation anyway though.
		//THE SOLUTION:
			//Define the variables as static, in the source file, WITHIN A PRIVATE MEMBER FUNCTION.
			//Thus the static variables themselves only have function scope.
			//Access is only available via the private member function, thus access is fully controlled.
			//They are shared amongst threads, so locks are necessary, but since they are private this class can handle it internally

		map<string, int>& Get_NumWritersByFileMap(void) const;
		map<string, size_t>& Get_FundamentalArraySizeMap(TTree* locTree) const;

		/************************************************************** MEMBER VARIABLES ************************************************************/

		TTree* dTree;
		string dFileName;
		string dTreeIndex_MajorBranchName;
		string dTreeIndex_MinorBranchName;

		/******************************************************** BRANCH MEMORY AND TYPE MAPPING ****************************************************/

		//These objects are kept here because they must be kept somewhere: 
			//branches addresses of pointers, so pointers must reside somewhere permanent
			//However, for fundamental objects/arrays: memory stored in the branches themselves: don't need to hold onto them
		map<string, TClonesArray*> dMemoryMap_ClonesArray;
		map<string, TObject*> dMemoryMap_TObject;
};

/******************************************************************** GET POINTERS ********************************************************************/

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

/******************************************************************* CREATE BRANCHES ******************************************************************/

template <typename DType> inline typename enable_if<std::is_base_of<TObject, DType>::value, void>::type
DTreeInterface::Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName)
{
	if(locArraySize == 0)
		Create_Branch_TObject<DType>(locBranchName);
	else
		Create_Branch_ClonesArray<DType>(locBranchName, locArraySize);
}

template <typename DType> inline typename enable_if<!std::is_base_of<TObject, DType>::value, void>::type
DTreeInterface::Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName)
{
	if(locArraySize == 0)
		Create_Branch_Fundamental<DType>(locBranchName);
	else
		Create_Branch_FundamentalArray<DType>(locBranchName, locArraySizeName, locArraySize);
}

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
}

template <typename DType> inline void DTreeInterface::Create_Branch_ClonesArray(string locBranchName, unsigned int locInitialSize)
{
	if(dTree->GetBranch(locBranchName.c_str()) != NULL)
		return; //already created

	dMemoryMap_ClonesArray[locBranchName] = new TClonesArray(DType::Class()->GetName(), locInitialSize);
	dTree->Branch(locBranchName.c_str(), &(dMemoryMap_ClonesArray[locBranchName]), 32000, 0); //0: don't split
}

/************************************************************** TEMPLATE SPECIALIZATIONS **************************************************************/

//template<> struct DTreeInterface::DROOTTypeString<const char*> { static const char* GetTypeString() {return "C";} };
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

/******************************************************************* MISCELLANEOUS ********************************************************************/

inline void DTreeInterface::Set_TreeIndexBranchNames(string locTreeIndex_MajorBranchName, string locTreeIndex_MinorBranchName)
{
	dTreeIndex_MajorBranchName = locTreeIndex_MajorBranchName;
	dTreeIndex_MinorBranchName = locTreeIndex_MinorBranchName;
}

//INCREASE ARRAY SIZE
template <typename DType> inline void DTreeInterface::Increase_ArraySize(string locBranchName, int locNewArraySize)
{
	//create a new, larger array if the current one is too small
		//DOES NOT copy the old results!  In other words, only call BETWEEN entries, not DURING an entry
	DType* locOldBranchAddress = Get_Pointer_Fundamental<DType>(locBranchName);
	dTree->SetBranchAddress(locBranchName.c_str(), new DType[locNewArraySize]);
	delete[] locOldBranchAddress;
}

template <typename DType> inline void DTreeInterface::Fill_TObject(string locBranchName, DType& locObject, bool locIsArrayFlag, size_t locArrayIndex)
{
	if(locIsArrayFlag)
	{
		TClonesArray* locClonesArray = Get_Pointer_TClonesArray(locBranchName);
		*(DType*)locClonesArray->ConstructedAt(locArrayIndex) = locObject;
	}
	else
		*Get_Pointer_TObject<DType>(locBranchName) = locObject;
}

#endif //DTreeInterface_h
