#ifndef DTreeInterfaceObjects_h
#define DTreeInterfaceObjects_h

#include <typeindex>
#include <map>
#include <string>

using namespace std;

class DTreeInterface;

/***************************************************************** DTreeTypeChecker *******************************************************************/

class DTreeTypeChecker
{
	template <typename DType> static void Is_Supported(void)
	{
		static_assert(std::is_same<Char_t, DType>::value || std::is_same<UChar_t, DType>::value ||
			std::is_same<Short_t, DType>::value || std::is_same<UShort_t, DType>::value || std::is_same<Int_t, DType>::value ||
			std::is_same<UInt_t, DType>::value || std::is_same<Float_t, DType>::value || std::is_same<Double_t, DType>::value ||
			std::is_same<Long64_t, DType>::value || std::is_same<ULong64_t, DType>::value || std::is_same<Bool_t, DType>::value ||
			std::is_same<TVector2, DType>::value || std::is_same<TVector3, DType>::value || std::is_same<TLorentzVector, DType>::value);
	}

	template <typename DType> static void Is_Fundamental(void)
	{
		static_assert(std::is_same<Char_t, DType>::value || std::is_same<UChar_t, DType>::value ||
			std::is_same<Short_t, DType>::value || std::is_same<UShort_t, DType>::value || std::is_same<Int_t, DType>::value ||
			std::is_same<UInt_t, DType>::value || std::is_same<Float_t, DType>::value || std::is_same<Double_t, DType>::value ||
			std::is_same<Long64_t, DType>::value || std::is_same<ULong64_t, DType>::value || std::is_same<Bool_t, DType>::value);
	}

	template <typename DType> static void Is_TObject(void)
	{
		static_assert(std::is_same<TVector2, DType>::value || std::is_same<TVector3, DType>::value || std::is_same<TLorentzVector, DType>::value);
	}
};

/**************************************************************** DTreeBranchRegister *****************************************************************/

class DTreeBranchRegister
{
	friend class DTreeInterface;

	public:
		template <typename DType> void Register_Branch_Single(string locBranchName);
		template <typename DType> void Register_Branch_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize = 10);
		template <typename DType> void Register_Branch_ClonesArray(string locBranchName, size_t locInitialArraySize = 10);

	private:
		map<string, type_index> dBranchTypeMap;
		map<string, size_t> dInitialArraySizeMap;
		map<string, string> dArraySizeNameMap; //for fundamental
};

template <typename DType> inline void DTreeBranchRegister::Register_Branch_Single(string locBranchName)
{
	DTreeTypeChecker::Is_Supported<DType>();
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
}

template <typename DType> inline void DTreeBranchRegister::Register_Branch_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize)
{
	DTreeTypeChecker::Is_Fundamental<DType>();
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
	dInitialArraySizeMap[locBranchName] = locInitialArraySize;
	dArraySizeNameMap[locBranchName] = locArraySizeName;
}

template <typename DType> inline void DTreeBranchRegister::Register_Branch_ClonesArray(string locBranchName, size_t locInitialArraySize)
{
	DTreeTypeChecker::Is_TObject<DType>();
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
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
		~DTreeFillData(void);
		template <typename DType> void Fill_Single(string locBranchName, const DType& locData);
		template <typename DType> void Fill_Array(string locBranchName, const DType& locData, unsigned int locArrayIndex);

	private:
		template <typename DType> void Fill_FundamentalArray(string locBranchName, const DType& locValue, unsigned int locArrayIndex);
		template <typename DType> void Fill_ClonesArray(string locBranchName, const DType& locObject, unsigned int locArrayIndex);

		void Delete(type_index locTypeIndex, void* locVoidPointer, bool locIsArrayFlag);
		template <typename DType> inline void Delete(DType* locObject, bool locIsArrayFlag);

		map<string, pair<type_index, void*> > dFillData;

		map<string, size_t> dArraySizeMap;
		map<string, size_t> dArrayNumFilledMap; //can be lest than the size (e.g. size of 10, but only first 3 filled)
};

/*********************************************************** DTreeFillData: FILL BRANCHES *************************************************************/

template <typename DType> inline void DTreeFillData::Fill_Single(string locBranchName, const DType& locData)
{
	DTreeTypeChecker::Is_Supported<DType>();
	type_index locTypeIndex(typeid(DType));

	map<string, pair<type_index, void*> >::iterator locIterator = dFillData.find(locBranchName);
	if(locIterator == dFillData.end())
	{
		dFillData[locBranchName] = pair<type_index, void*>(locTypeIndex, new DType(locData));
		return;
	}

	if(locTypeIndex != locIterator->second.first)
	{
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
		return;
	}

	void* locVoid = locIterator->second.second;
	*(static_cast<DType*>(locVoid)) = locData;
}

template <typename DType> inline void DTreeFillData::Fill_Array(string locBranchName, const DType& locData, unsigned int locArrayIndex)
{
	DTreeTypeChecker::Is_Supported<DType>();
	if(std::is_base_of<TObject, DType>::value)
		Fill_FundamentalArray<DType>(locBranchName, locData, locArrayIndex);
	else
		Fill_ClonesArray<DType>(locBranchName, locData, locArrayIndex);
}

template <typename DType> inline void DTreeFillData::Fill_FundamentalArray(string locBranchName, const DType& locData, size_t locArrayIndex)
{
	type_index locTypeIndex(typeid(DType));

	map<string, pair<type_index, void*> >::iterator locIterator = dFillData.find(locBranchName);
	void* locVoid = NULL;

	if(locIterator == dFillData.end())
	{
		locVoid = new DType[10];
		dFillData[locBranchName] = pair<type_index, void*>(locTypeIndex, locVoid);
		dArraySizeMap[locBranchName] = 10;
	}
	else if(locTypeIndex != locIterator->second.first)
	{
		cout << "WARNING: CANNOT FILL: IS WRONG TYPE FOR BRANCH " << locBranchName << endl;
		return;
	}
	else
		locVoid = locIterator->second.second;

	DType* locArray = static_cast<DType*>(locVoid);

	//create a new, larger array if the current one is too small
	unsigned int locCurrentArraySize = dArraySizeMap[locBranchName];
	if(locArrayIndex >= locCurrentArraySize)
	{
		DType* locOldArray = locArray;
		locArray = new DType[locArrayIndex + 1];
		dFillData[locBranchName].second = locArray;

		//copy the old contents into the new array
		for(unsigned int loc_i = 0; loc_i < locCurrentArraySize; ++loc_i)
			locArray[loc_i] = locOldArray[loc_i];

		delete[] locOldArray;
		dArraySizeMap[locBranchName] = locArrayIndex + 1;
	}

	//set the data
	locArray[locArrayIndex] = locData;
}

template <typename DType> inline void DTreeFillData::Fill_ClonesArray(string locBranchName, const DType& locObject, unsigned int locArrayIndex)
{
	TClonesArray* locClonesArray = Get_Pointer_TClonesArray(locBranchName);
	*(DType*)locClonesArray->ConstructedAt(locArrayIndex) = locObject;
}

/************************************************************* DTreeFillData: DESTRUCTOR **************************************************************/

inline DTreeFillData::~DTreeFillData(void)
{
	//delete all memory (void*'s)
	//loop over branches
	for(auto locBranchIterator : dFillData)
	{
		string locBranchName = locBranchIterator->first;
		type_index locTypeIndex = locBranchIterator->second.first;
		void* locVoidPointer = locBranchIterator->second.second;
		bool locIsArrayFlag = (dArrayNumFilledMap.find(locBranchName) != locTreeFillData->dArrayNumFilledMap.end());

		Delete(locTypeIndex, locVoidPointer, locIsArrayFlag);
	}
}

inline void DTreeFillData::Delete(type_index locTypeIndex, void* locVoidPointer, bool locIsArrayFlag)
{
	//Fundamental types
	if(locTypeIndex == type_index(typeid(Char_t)))
		Delete(static_cast<Char_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(UChar_t)))
		Delete(static_cast<UChar_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Short_t)))
		Delete(static_cast<Short_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(UShort_t)))
		Delete(static_cast<UShort_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Int_t)))
		Delete(static_cast<Int_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(UInt_t)))
		Delete(static_cast<UInt_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Float_t)))
		Delete(static_cast<Float_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Double_t)))
		Delete(static_cast<Double_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Long64_t)))
		Delete(static_cast<Long64_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(ULong64_t)))
		Delete(static_cast<ULong64_t>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(Bool_t)))
		Delete(static_cast<Bool_t>(locVoidPointer), locIsArrayFlag);

	//TObject
	else if(locTypeIndex == type_index(typeid(TVector3)))
		Delete(static_cast<TVector3>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(TVector2)))
		Delete(static_cast<TVector2>(locVoidPointer), locIsArrayFlag);
	else if(locTypeIndex == type_index(typeid(TLorentzVector)))
		Delete(static_cast<TLorentzVector>(locVoidPointer), locIsArrayFlag);
}

template <typename DType> inline void DTreeFillData::Delete(DType* locObject, bool locIsArrayFlag)
{
	if(locIsArrayFlag)
		delete[] locObject;
	else
		delete locObject;		
}

#endif //DTreeInterfaceObjects
