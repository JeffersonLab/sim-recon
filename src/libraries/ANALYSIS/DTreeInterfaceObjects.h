#ifndef DTreeInterfaceObjects_h
#define DTreeInterfaceObjects_h

#include <typeindex>
#include <typeinfo>
#include <map>
#include <string>
#include <deque>
#include <vector>

#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

class DTreeInterface;

/***************************************************************** DTreeTypeChecker *******************************************************************/

class DTreeTypeChecker
{
	public:
		template <typename DType> static void Is_Supported(void)
		{
			static_assert(std::is_same<Char_t, DType>::value || std::is_same<UChar_t, DType>::value ||
				std::is_same<Short_t, DType>::value || std::is_same<UShort_t, DType>::value || std::is_same<Int_t, DType>::value ||
				std::is_same<UInt_t, DType>::value || std::is_same<Float_t, DType>::value || std::is_same<Double_t, DType>::value ||
				std::is_same<Long64_t, DType>::value || std::is_same<ULong64_t, DType>::value || std::is_same<Bool_t, DType>::value ||
				std::is_same<TVector2, DType>::value || std::is_same<TVector3, DType>::value || std::is_same<TLorentzVector, DType>::value,
				"DTreeTypeChecker ERROR: TYPE IS NOT SUPPORTED.");
		}

		template <typename DType> static void Is_Fundamental(void)
		{
			static_assert(std::is_same<Char_t, DType>::value || std::is_same<UChar_t, DType>::value ||
				std::is_same<Short_t, DType>::value || std::is_same<UShort_t, DType>::value || std::is_same<Int_t, DType>::value ||
				std::is_same<UInt_t, DType>::value || std::is_same<Float_t, DType>::value || std::is_same<Double_t, DType>::value ||
				std::is_same<Long64_t, DType>::value || std::is_same<ULong64_t, DType>::value || std::is_same<Bool_t, DType>::value,
				"DTreeTypeChecker ERROR: TYPE IS NOT A SUPPORTED FUNDAMENTAL TYPE.");
		}

		template <typename DType> static void Is_TObject(void)
		{
			static_assert(std::is_same<TVector2, DType>::value || std::is_same<TVector3, DType>::value || std::is_same<TLorentzVector, DType>::value,
					"DTreeTypeChecker ERROR: TYPE IS NOT A SUPPORTED TOBJECT TYPE.");
		}
};

/**************************************************************** DTreeBranchRegister *****************************************************************/

class DTreeBranchRegister
{
	friend class DTreeInterface;

	public:
		DTreeBranchRegister(void) : dUserInfo(new TList()){}
		~DTreeBranchRegister(void)
		{
			//If multiple threads, objects in the list may be created by multiple threads on the heap, each with the same name
			//The TList destructor tries to do something I don't quite understand ...
			//It looks like it's trying to delete objects on the heap, which it shouldn't try to do by default
			//Anyway, it looks like things are being double-freed somehow
			//Instead: Avoid this by removing the list entries
			while(dUserInfo->GetEntries() > 0)
				dUserInfo->RemoveLast();
			//delete dUserInfo; //in fact, deleting this STILL seems to be causing issues.  Just don't delete it. 
		}

		TList* Get_UserInfo(void) const{return dUserInfo;}

		template <typename DType> void Register_Single(string locBranchName);
		template <typename DType> void Register_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize = 10);
		template <typename DType> void Register_ClonesArray(string locBranchName, size_t locInitialArraySize = 10);

	private:
		TList* dUserInfo;
		vector<string> dBranchNames; //keep the order in which they were added
		map<string, type_index> dBranchTypeMap;
		map<string, size_t> dInitialArraySizeMap;
		map<string, string> dArraySizeNameMap; //for fundamental
};

template <typename DType> inline void DTreeBranchRegister::Register_Single(string locBranchName)
{
	DTreeTypeChecker::Is_Supported<DType>();
	dBranchNames.push_back(locBranchName);
	dBranchTypeMap.insert(pair<string, type_index>(locBranchName, type_index(typeid(DType))));
}

template <typename DType> inline void DTreeBranchRegister::Register_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize)
{
	DTreeTypeChecker::Is_Fundamental<DType>();
	dBranchNames.push_back(locBranchName);
	dBranchTypeMap.insert(pair<string, type_index>(locBranchName, type_index(typeid(DType))));
	dInitialArraySizeMap[locBranchName] = locInitialArraySize;
	dArraySizeNameMap[locBranchName] = locArraySizeName;
}

template <typename DType> inline void DTreeBranchRegister::Register_ClonesArray(string locBranchName, size_t locInitialArraySize)
{
	DTreeTypeChecker::Is_TObject<DType>();
	dBranchNames.push_back(locBranchName);
	dBranchTypeMap.insert(pair<string, type_index>(locBranchName, type_index(typeid(DType))));
	dInitialArraySizeMap[locBranchName] = locInitialArraySize;
}

/******************************************************************* DTreeFillData ********************************************************************/

//Want to abstract the fill type so we can hold them in a container without dynamically allocating void*'s
class DFillBaseClass
{
	public:
		virtual ~DFillBaseClass(){};
		virtual void* Get(size_t locArrayIndex) = 0;
		virtual void Check_Capacity(void) = 0;
};

template <typename DType>
class DFillClass : public DFillBaseClass
{
	public:
		deque<DType> dFillData; //must use deque because vector<bool> won't compile!!!

		~DFillClass(){};

		void* Get(size_t locArrayIndex){return static_cast<void*>(&(dFillData[locArrayIndex]));}
		void Check_Capacity(void);

	private:
		size_t dMaxFillVectorSize = 1000; //if exceeds this, will drop down on next event
};

template <typename DType> inline void DFillClass<DType>::Check_Capacity(void)
{
	if(dFillData.size() <= dMaxFillVectorSize)
		return;
	dFillData.resize(dMaxFillVectorSize);
}

//Need one per thread:
	//If this is created within the scope of a single object that is shared amongst all threads (e.g. plugin processor): static thread-local variable
		//Data stored as void*: Requires new on creation and delete on destruction: Try to re-use object
	//If this is created within the scope of an object that is unique for each thread: Class member
class DTreeFillData
{
	friend class DTreeInterface;

	public:
		~DTreeFillData(void);
		template <typename DType> void Fill_Single(string locBranchName, const DType& locData);
		template <typename DType> void Fill_Array(string locBranchName, const DType& locData, size_t locArrayIndex);

	private:
		map<string, pair<type_index, DFillBaseClass*> > dFillData;
		map<string, int> dArrayLargestIndexFilledMap; //can be less than the size //reset by DTreeInterface after fill
};

/*********************************************************** DTreeFillData: FILL BRANCHES *************************************************************/

template <typename DType> inline void DTreeFillData::Fill_Single(string locBranchName, const DType& locData)
{
	DTreeTypeChecker::Is_Supported<DType>();
	type_index locTypeIndex(typeid(DType));

	auto locIterator = dFillData.find(locBranchName);
	if(locIterator == dFillData.end())
	{
		//create new object, register in map
		auto locFillClass = new DFillClass<DType>();
		locFillClass->dFillData.push_back(locData);
		dFillData.emplace(locBranchName, std::make_pair(locTypeIndex, static_cast<DFillBaseClass*>(locFillClass)));
	}
	else if(locTypeIndex != locIterator->second.first)
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
	else
	{
		auto locFillClass = static_cast<DFillClass<DType>*>(locIterator->second.second);
		locFillClass->dFillData[0] = locData;
	}
}

template <typename DType> inline void DTreeFillData::Fill_Array(string locBranchName, const DType& locData, size_t locArrayIndex)
{
	DTreeTypeChecker::Is_Supported<DType>();
	type_index locTypeIndex(typeid(DType));

	auto locIterator = dFillData.find(locBranchName);
	if(locIterator == dFillData.end())
	{
		//create new object, register in map
		auto locFillClass = new DFillClass<DType>();
		dFillData.emplace(locBranchName, std::make_pair(locTypeIndex, static_cast<DFillBaseClass*>(locFillClass)));

		//fill
		locFillClass->dFillData.resize(locArrayIndex + 1);
		locFillClass->dFillData[locArrayIndex] = locData;
		dArrayLargestIndexFilledMap[locBranchName] = locArrayIndex;
	}
	else if(locTypeIndex != locIterator->second.first)
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
	else
	{
		auto locFillClass = static_cast<DFillClass<DType>*>(locIterator->second.second);

		//resize if needed & fill
		if(locArrayIndex >= locFillClass->dFillData.size())
			locFillClass->dFillData.resize(locArrayIndex + 1);
		locFillClass->dFillData[locArrayIndex] = locData;

		//register largest index filled
		auto& locLargestIndexFilled = dArrayLargestIndexFilledMap[locBranchName];
		if(int(locArrayIndex) > locLargestIndexFilled)
			locLargestIndexFilled = locArrayIndex;
	}
}
/*
//Enable this version if type inherits from TObject //void: is return type
template <typename DType> typename enable_if<std::is_base_of<TObject, DType>::value, void>::type
		Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName);

//Enable this version if type does NOT inherit from TObject //void: is return type
template <typename DType> typename enable_if<!std::is_base_of<TObject, DType>::value, void>::type
		Create_Branch(string locBranchName, size_t locArraySize, string locArraySizeName);
*/
/************************************************************* DTreeFillData: DESTRUCTOR **************************************************************/

inline DTreeFillData::~DTreeFillData(void)
{
	//delete all memory (void*'s)
	//loop over branches
	for(auto& locBranchPair : dFillData)
		delete locBranchPair.second.second;
}

#endif //DTreeInterfaceObjects
