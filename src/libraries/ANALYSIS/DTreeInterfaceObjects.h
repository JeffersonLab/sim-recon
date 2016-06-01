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
		~DTreeBranchRegister(void){delete dUserInfo;}

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

//Need one per thread:
	//If this is created within the scope of a single object that is shared amongst all threads (e.g. plugin processor): static thread-local variable
		//Data stored as void*: Requires new on creation and delete on destruction: Try to re-use object
	//If this is created within the scope of an object that is unique for each thread: Class member
class DTreeFillData
{
	friend class DTreeInterface;

	public:
		DTreeFillData(void) : dFillData(new map<string, pair<type_index, deque<void*>* > >()), dArrayLargestIndexFilledMap(new map<string, size_t>()){}
		~DTreeFillData(void);
		template <typename DType> void Fill_Single(string locBranchName, const DType& locData);
		template <typename DType> void Fill_Array(string locBranchName, const DType& locData, unsigned int locArrayIndex);

	private:
		void Delete(type_index& locTypeIndex, deque<void*>& locVoidDeque);

		map<string, pair<type_index, deque<void*>* > >* dFillData;
		map<string, size_t>* dArrayLargestIndexFilledMap; //can be less than the size //reset by DTreeInterface after fill
};

/*********************************************************** DTreeFillData: FILL BRANCHES *************************************************************/

template <typename DType> inline void DTreeFillData::Fill_Single(string locBranchName, const DType& locData)
{
	DTreeTypeChecker::Is_Supported<DType>();
	type_index locTypeIndex(typeid(DType));

	auto locIterator = dFillData->find(locBranchName);
	if(locIterator == dFillData->end())
	{
		void* locVoidData = static_cast<void*>(new DType(locData));
		deque<void*>* locVoidDeque = new deque<void*>(1, locVoidData);
		pair<type_index, deque<void*>* > locTypePair(locTypeIndex, locVoidDeque);
		pair<string, pair<type_index, deque<void*>* > > locMapPair(locBranchName, locTypePair);
		dFillData->insert(locMapPair);
		return;
	}

	if(locTypeIndex != locIterator->second.first)
	{
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
		return;
	}

	deque<void*>& locVoidDeque = *(locIterator->second.second);
	*(static_cast<DType*>(locVoidDeque[0])) = locData;
}

template <typename DType> inline void DTreeFillData::Fill_Array(string locBranchName, const DType& locData, unsigned int locArrayIndex)
{
	DTreeTypeChecker::Is_Supported<DType>();
	type_index locTypeIndex(typeid(DType));

	auto locIterator = dFillData->find(locBranchName);
	if(locIterator == dFillData->end())
	{
		void* locVoidData = static_cast<void*>(new DType(locData));
		deque<void*>* locVoidDeque = new deque<void*>(locArrayIndex + 1, nullptr);
		(*locVoidDeque)[locArrayIndex] = locVoidData;
		pair<type_index, deque<void*>* > locTypePair(locTypeIndex, locVoidDeque);
		pair<string, pair<type_index, deque<void*>* > > locMapPair(locBranchName, locTypePair);
		dFillData->insert(locMapPair);
		(*dArrayLargestIndexFilledMap)[locBranchName] = locArrayIndex;
		return;
	}
	else if(locTypeIndex != locIterator->second.first)
	{
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
		return;
	}

	deque<void*>& locVoidDeque = *(locIterator->second.second);

	//expand deque if needed
	for(size_t loc_i = locVoidDeque.size(); loc_i <= locArrayIndex; ++loc_i)
		locVoidDeque.push_back(static_cast<void*>(new DType));

	//set the data
	*(static_cast<DType*>(locVoidDeque[locArrayIndex])) = locData;

	//register largest index filled
	auto& locLargestIndexFilled = (*dArrayLargestIndexFilledMap)[locBranchName];
	if(locArrayIndex > locLargestIndexFilled)
		locLargestIndexFilled = locArrayIndex;
}

/************************************************************* DTreeFillData: DESTRUCTOR **************************************************************/

inline DTreeFillData::~DTreeFillData(void)
{
	//delete all memory (void*'s)
	//loop over branches
	for(auto locBranchIterator : *dFillData)
	{
		string locBranchName = locBranchIterator.first;
		type_index& locTypeIndex = locBranchIterator.second.first;
		deque<void*>* locVoidDeque = locBranchIterator.second.second;
		Delete(locTypeIndex, *locVoidDeque);
		delete locVoidDeque;
	}

	delete dFillData;
	delete dArrayLargestIndexFilledMap;
}

inline void DTreeFillData::Delete(type_index& locTypeIndex, deque<void*>& locVoidDeque)
{
	for(size_t loc_i = 0; loc_i < locVoidDeque.size(); ++loc_i)
	{
		//Fundamental types
		if(locTypeIndex == type_index(typeid(Char_t)))
			delete (static_cast<Char_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(UChar_t)))
			delete (static_cast<UChar_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Short_t)))
			delete (static_cast<Short_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(UShort_t)))
			delete (static_cast<UShort_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Int_t)))
			delete (static_cast<Int_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(UInt_t)))
			delete (static_cast<UInt_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Float_t)))
			delete (static_cast<Float_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Double_t)))
			delete (static_cast<Double_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Long64_t)))
			delete (static_cast<Long64_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(ULong64_t)))
			delete (static_cast<ULong64_t*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(Bool_t)))
			delete (static_cast<Bool_t*>(locVoidDeque[loc_i]));

		//TObject
		else if(locTypeIndex == type_index(typeid(TVector3)))
			delete (static_cast<TVector3*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(TVector2)))
			delete (static_cast<TVector2*>(locVoidDeque[loc_i]));
		else if(locTypeIndex == type_index(typeid(TLorentzVector)))
			delete (static_cast<TLorentzVector*>(locVoidDeque[loc_i]));
	}
}

#endif //DTreeInterfaceObjects
