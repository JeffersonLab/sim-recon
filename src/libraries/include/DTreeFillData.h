#include <typeindex>
#include <map>
#include <string>

using namespace std;

class DTreeInterface;

class DTreeBranchRegister
{
	friend class DTreeInterface;

	public:
		template <typename DType> void Register_Branch(string locBranchName);
		template <typename DType> void Register_Branch_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize = 10);
		template <typename DType> void Register_Branch_ClonesArray(string locBranchName, size_t locInitialArraySize = 10);

	private:
		map<string, type_index> dBranchTypeMap;
		map<string, size_t> dInitialArraySizeMap;
		map<string, string> dArraySizeNameMap; //for fundamental
};

template <typename DType> inline void DTreeBranchRegister::Register_Branch(string locBranchName)
{
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
}

template <typename DType> inline void DTreeBranchRegister::Register_Branch_FundamentalArray(string locBranchName, string locArraySizeName, size_t locInitialArraySize)
{
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
	dInitialArraySizeMap[locBranchName] = locInitialArraySize;
	dArraySizeNameMap[locBranchName] = locArraySizeName;
}

template <typename DType> inline void DTreeBranchRegister::Register_Branch_ClonesArray(string locBranchName, size_t locInitialArraySize)
{
	dBranchTypeMap[locBranchName] = type_index(typeid(DType));
	dInitialArraySizeMap[locBranchName] = locInitialArraySize;
}

//Need one per thread
	//In plugin: Local variable in evnt function (ugh)
		//Ugh: need to internally hold void*'s: have to call "new" every time
	//In action: Class member (doesn't work if put directly in plugin (ugh))
class DTreeFillData
{
	friend class DTreeInterface;

	public:
		template <typename DType> void Fill_Data(string locBranchName, DType locValue);
		template <typename DType> void Fill_Data(string locBranchName, DType locValue, unsigned int locArrayIndex);
//		template <typename DType> void Fill_TObject(string locBranchName, DType& locObject, unsigned int locArrayIndex);
//		template <typename DType> void Fill_TObject(string locBranchName, DType& locObject);

	private:
		map<string, pair<type_index, void*> > dFillData;

		map<string, size_t> dArraySizeMap;
		map<string, size_t> dArrayNumFilledMap; //can be lest than the size (e.g. size of 10, but only first 3 filled)
};

/******************************************************************* FILL BRANCHES ********************************************************************/

//FILL BRANCHES
template <typename DType> inline void DTreeFillData::Fill_Data(string locBranchName, DType& locData)
{
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

template <typename DType> inline void DTreeFillData::Fill_Array(string locBranchName, DType& locData, size_t locArrayIndex)
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

template <typename DType> inline void DTreeFillData::Fill_TObject(string locBranchName, DType& locObject, unsigned int locArrayIndex)
{
	TClonesArray* locClonesArray = Get_Pointer_TClonesArray(locBranchName);
	*(DType*)locClonesArray->ConstructedAt(locArrayIndex) = locObject;
}

