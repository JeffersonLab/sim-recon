#ifndef _DEventWriterROOT_
#define _DEventWriterROOT_

#include <map>
#include <string>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"

#include <JANA/JApplication.h>
#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>

#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DTrackTimeBased.h>

#include <PID/DChargedTrack.h>
#include <PID/DBeamPhoton.h>
#include <PID/DChargedTrackHypothesis.h>
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DNeutralShower.h>
#include <PID/DEventRFBunch.h>
#include <PID/DMCReaction.h>

#include <ANALYSIS/DParticleCombo.h>
#include <ANALYSIS/DReaction.h>
#include <ANALYSIS/DAnalysisResults.h>
#include <ANALYSIS/DAnalysisUtilities.h>
#include <ANALYSIS/DMCThrownMatching.h>
#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DEventWriterROOT : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterROOT);

		DEventWriterROOT(JEventLoop* locEventLoop);
		virtual ~DEventWriterROOT(void);

		void Create_DataTrees(JEventLoop* locEventLoop) const;
		void Create_ThrownTree(string locOutputFileName) const;

		void Fill_DataTrees(JEventLoop* locEventLoop, string locDReactionTag) const; //fills all from this factory tag
		void Fill_DataTree(JEventLoop* locEventLoop, const DReaction* locReaction, deque<const DParticleCombo*>& locParticleCombos) const;
		void Fill_ThrownTree(JEventLoop* locEventLoop) const;

	protected:

		//CUSTOM FUNCTIONS: //Inherit from this class and write custom code in these functions
			//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier: You don't need to manage the branch memory.
			//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
			//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
				//Note that the JEventLoop is unavailable.  This is to prevent calls to other factories that may cause deadlock.
			//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.
		virtual void Create_CustomBranches_ThrownTree(TTree* locTree) const{};
		virtual void Fill_CustomBranches_ThrownTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const{};
		virtual void Create_CustomBranches_DataTree(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag) const{};
		virtual void Fill_CustomBranches_DataTree(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
				const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
				const vector<const DNeutralParticle*>& locNeutralParticles, const deque<const DParticleCombo*>& locParticleCombos) const{};

		//UTILITY FUNCTIONS
		string Convert_ToBranchName(string locInputName) const;
		string Build_BranchName(string locParticleBranchName, string locVariableName) const;
		ULong64_t Calc_ParticleMultiplexID(Particle_t locPID) const;
		void Get_DecayProductNames(const DReaction* locReaction, size_t locReactionStepIndex, TMap* locPositionToNameMap, TList*& locDecayProductNames, deque<size_t>& locSavedSteps) const;

		//BRANCH CREATION: //with the full branch name
		template <typename DType> void Create_Branch_Fundamental(TTree* locTree, string locBranchName) const;
		template <typename DType> void Create_Branch_NoSplitTObject(TTree* locTree, string locBranchName) const;
		template <typename DType> void Create_Branch_FundamentalArray(TTree* locTree, string locBranchName, string locArraySizeString, unsigned int locInitialSize) const;
		void Create_Branch_ClonesArray(TTree* locTree, string locBranchName, string locClassName, unsigned int locSize) const;

		//BRANCH CREATION: //with separate particle & variable names (from which the branch name is made)
		template <typename DType> void Create_Branch_Fundamental(TTree* locTree, string locParticleBranchName, string locVariableName) const;
		template <typename DType> void Create_Branch_NoSplitTObject(TTree* locTree, string locParticleBranchName, string locVariableName) const;
		template <typename DType> void Create_Branch_FundamentalArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locArraySizeString, unsigned int locInitialSize) const;
		void Create_Branch_ClonesArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locClassName, unsigned int locSize) const;

		//BRANCH FILLING: //with the full branch name
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue) const;
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue, unsigned int locArrayIndex) const;
		template <typename DType> void Fill_ClonesData(TTree* locTree, string locBranchName, DType& locObject, unsigned int locArrayIndex) const;
		template <typename DType> void Fill_TObjectData(TTree* locTree, string locBranchName, DType& locObject) const;

		//BRANCH FILLING: //with separate particle & variable names (from which the branch name is made)
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue) const;
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex) const;
		template <typename DType> void Fill_ClonesData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject, unsigned int locArrayIndex) const;
		template <typename DType> void Fill_TObjectData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject) const;

		const DAnalysisUtilities* dAnalysisUtilities;

	private:

		unsigned int dInitNumThrownArraySize;
		unsigned int dInitNumBeamArraySize;
		unsigned int dInitNumTrackArraySize;
		unsigned int dInitNumNeutralArraySize;
		unsigned int dInitNumComboArraySize;

		string dTrackSelectionTag;
		string dShowerSelectionTag;

		//DEFAULT ACTIONS LISTED SEPARATELY FROM CUSTOM (in case in derived class user does something bizarre)
		map<const DReaction*, DCutAction_ThrownTopology*> dCutActionMap_ThrownTopology;
		map<const DReaction*, DCutAction_TrueCombo*> dCutActionMap_TrueCombo;
		map<const DReaction*, DCutAction_BDTSignalCombo*> dCutActionMap_BDTSignalCombo;

		//add in future: let user execute custom actions (outside of lock): user adds and initializes actions in derived-writer constructor
		//map<const DReaction*, map<string, map<const DParticleCombo*, bool> > > dCustomActionResultsMap; //string is action name

		DEventWriterROOT(void){}; //don't allow default constructor

		/****************************************** STATIC-VARIABLE-ACCESSING PRIVATE MEMBER FUNCTIONS ******************************************/

		//Some variables needs to be shared amongst threads (e.g. the memory used for the branch variables)
			//They must be (indirectly) accessible to derived classes (custom branches in derived classes)
		//However, you cannot make them global/extern/static/static-member variables in the header file:
			//They would be in the header file, and the header file is included in the ANALYSIS library AND in each plugin that uses it
				//When a header file is included in a src file, it's contents are essentially copied directly into it
			//Thus there are two instances of each static variable: one in each translation unit (library)
			//Supposedly(?) they are linked together during runtime when loading, so there is (supposedly) no undefined behavior.
			//However, this causes a double free (double-deletion) when these libraries are closed at the end of the program, crashing it.
		//Thus the variables must be in a single source file that is compiled into a single library
		//However, you (somehow?) cannot make them global/extern variables in the cpp function
			//This also (somehow?) causes the double-free problem above for (at least) stl containers
			//It works for pointers-to-stl-containers and fundamental types, but I dunno why.
			//It's not good encapsulation anyway though.
		//THE SOLUTION:
			//Define the variables as static, in the source file, WITHIN A PRIVATE MEMBER FUNCTION.
			//Thus the static variables themselves only have function scope.
			//Access is only available via the private member function, thus access is fully controlled.
			//They are shared amongst threads, so locks are necessary, but since they are private this class can handle it internally

		int& Get_NumEventWriterThreads(void) const;

		//keep track of thrown tree & file name
		pair<string, TTree*>& Get_ThrownTreePair(void) const; //string is file name

		//keep track of the array sizes used for the branch memory
		map<TTree*, map<string, unsigned int> >& Get_FundamentalArraySizeMap(void) const;

		//keep track of the objects used for the branch memory
		map<TTree*, map<string, TClonesArray*> >& Get_ClonesArrayMap(void) const; //string is branch name
		map<TTree*, map<string, TObject*> >& Get_TObjectMap(void) const; //string is branch name

		//keep track of the output files
		map<string, TFile*>& Get_OutputROOTFileMap(void) const;

		//keep track of the ttrees
		map<string, TTree*>& Get_TTreeMap(void) const;

		/****************************************************************************************************************************************/

		void Get_Reactions(jana::JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const;

		//TREE CREATION:
		void Create_DataTree(const DReaction* locReaction, bool locIsMCDataFlag, double locTargetCenterZ) const;
		void Create_UserInfoMaps(TTree* locTree, const DReaction* locReaction, map<Particle_t, unsigned int>& locParticleNumberMap, double locTargetCenterZ) const;
		void Create_UserTargetInfo(TTree* locTree, Particle_t locTargetPID, double locTargetCenterZ) const;
		void Create_Branches_Thrown(TTree* locTree, bool locIsOnlyThrownFlag) const;

		//TREE CREATION: PARTICLE INFO
		void Create_Branches_ThrownParticles(TTree* locTree, bool locIsOnlyThrownFlag) const;
		void Create_Branches_Beam(TTree* locTree, bool locIsMCDataFlag) const;
		void Create_Branches_NeutralHypotheses(TTree* locTree, bool locIsMCDataFlag) const;
		void Create_Branches_ChargedHypotheses(TTree* locTree, bool locIsMCDataFlag) const;

		//TREE CREATION: COMBO INFO
		void Create_Branches_Combo(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag, const map<Particle_t, unsigned int>& locParticleNumberMap) const;
		void Create_Branches_BeamComboParticle(TTree* locTree, Particle_t locBeamPID, DKinFitType locKinFitType) const;
		void Create_Branches_ComboTrack(TTree* locTree, string locParticleBranchName, DKinFitType locKinFitType) const;
		void Create_Branches_ComboNeutral(TTree* locTree, string locParticleBranchName, DKinFitType locKinFitType) const;

		//TREE FILLING: THROWN INFO
		void Compute_ThrownPIDInfo(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying,
				ULong64_t& locNumPIDThrown_FinalState, ULong64_t& locPIDThrown_Decaying) const;
		void Group_ThrownParticles(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying,
				vector<const DMCThrown*>& locMCThrownsToSave, map<const DMCThrown*, unsigned int>& locThrownIndexMap) const;
		void Fill_ThrownInfo(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying) const;
		void Fill_ThrownInfo(TTree* locTree, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying,
				const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const;
		void Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, const DMCThrown* locMCThrown, const map<const DMCThrown*, unsigned int>& locThrownIndexMap,
				const DMCThrownMatching* locMCThrownMatching, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const;

		//TREE FILLING: INDEPENDENT PARTICLES
		void Fill_BeamData(TTree* locTree, unsigned int locArrayIndex, const DBeamPhoton* locBeamPhoton, const DVertex* locVertex, const DMCThrownMatching* locMCThrownMatching) const;
		void Fill_ChargedHypo(TTree* locTree, unsigned int locArrayIndex, const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const;
		void Fill_NeutralHypo(TTree* locTree, unsigned int locArrayIndex, const DNeutralParticleHypothesis* locPhotonHypothesis, const DMCThrownMatching* locMCThrownMatching,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const;

		//TREE FILLING: COMBO
		void Fill_ComboData(TTree* locTree, const DParticleCombo* locParticleCombo, unsigned int locComboIndex, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const;
		void Fill_ComboStepData(TTree* locTree, const DParticleCombo* locParticleCombo, unsigned int locStepIndex, unsigned int locComboIndex,
				DKinFitType locKinFitType, const map<string, map<oid_t, int> >& locObjectToArrayIndexMap) const;

		//TREE FILLING: COMBO PARTICLES
		void Fill_ComboBeamData(TTree* locTree, unsigned int locComboIndex, const DBeamPhoton* locBeamPhoton, unsigned int locBeamIndex, DKinFitType locKinFitType) const;
		void Fill_ComboChargedData(TTree* locTree, unsigned int locComboIndex, string locParticleBranchName, const DChargedTrackHypothesis* locMeasuredChargedHypo,
				const DChargedTrackHypothesis* locChargedHypo, unsigned int locChargedIndex, DKinFitType locKinFitType) const;
		void Fill_ComboNeutralData(TTree* locTree, unsigned int locComboIndex, string locParticleBranchName, const DNeutralParticleHypothesis* locMeasuredNeutralHypo,
				const DNeutralParticleHypothesis* locNeutralHypo, unsigned int locNeutralIndex, DKinFitType locKinFitType) const;

		//For ROOT type string for fundamental data variables
			//Defined in https://root.cern.ch/root/htmldoc/TTree.html
		// Main template class
		template<typename DType> struct DROOTTypeString { static const char* GetTypeString() {return "";} };
};

// Template specializations
template<> struct DEventWriterROOT::DROOTTypeString<const char*> { static const char* GetTypeString() {return "C";} };
template<> struct DEventWriterROOT::DROOTTypeString<Char_t> { static const char* GetTypeString() {return "B";} };
template<> struct DEventWriterROOT::DROOTTypeString<UChar_t> { static const char* GetTypeString() {return "b";} };
template<> struct DEventWriterROOT::DROOTTypeString<Short_t> { static const char* GetTypeString() {return "S";} };
template<> struct DEventWriterROOT::DROOTTypeString<UShort_t> { static const char* GetTypeString() {return "s";} };
template<> struct DEventWriterROOT::DROOTTypeString<Int_t> { static const char* GetTypeString() {return "I";} };
template<> struct DEventWriterROOT::DROOTTypeString<UInt_t> { static const char* GetTypeString() {return "i";} };
template<> struct DEventWriterROOT::DROOTTypeString<Float_t> { static const char* GetTypeString() {return "F";} };
template<> struct DEventWriterROOT::DROOTTypeString<Double_t> { static const char* GetTypeString() {return "D";} };
template<> struct DEventWriterROOT::DROOTTypeString<Long64_t> { static const char* GetTypeString() {return "L";} };
template<> struct DEventWriterROOT::DROOTTypeString<ULong64_t> { static const char* GetTypeString() {return "l";} };
template<> struct DEventWriterROOT::DROOTTypeString<Bool_t> { static const char* GetTypeString() {return "O";} };

//BRANCH CREATION: //with separate particle & variable names (from which the branch name is made)
template <typename DType> void DEventWriterROOT::Create_Branch_Fundamental(TTree* locTree, string locParticleBranchName, string locVariableName) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Create_Branch_Fundamental<DType>(locTree, locBranchName);
}

template <typename DType> void DEventWriterROOT::Create_Branch_NoSplitTObject(TTree* locTree, string locParticleBranchName, string locVariableName) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Create_Branch_NoSplitTObject<DType>(locTree, locBranchName);
}

template <typename DType> void DEventWriterROOT::Create_Branch_FundamentalArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locArraySizeString, unsigned int locInitialSize) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Create_Branch_FundamentalArray<DType>(locTree, locBranchName, locArraySizeString, locInitialSize);
}

inline void DEventWriterROOT::Create_Branch_ClonesArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locClassName, unsigned int locSize) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Create_Branch_ClonesArray(locTree, locBranchName, locClassName, locSize);
}

//BRANCH CREATION: //with the full branch name
template <typename DType> void DEventWriterROOT::Create_Branch_Fundamental(TTree* locTree, string locBranchName) const
{
	string locTypeString = DROOTTypeString<DType>::GetTypeString();
	string locTypeName = locBranchName + string("/") + locTypeString;
	locTree->Branch(locBranchName.c_str(), new DType(), locTypeName.c_str());
}

template <typename DType> void DEventWriterROOT::Create_Branch_NoSplitTObject(TTree* locTree, string locBranchName) const
{
	Get_TObjectMap()[locTree].insert(pair<string, TObject*>(locBranchName, (TObject*)(new DType())));
	locTree->Branch<DType>(locBranchName.c_str(), (DType**)&(Get_TObjectMap()[locTree][locBranchName]), 32000, 0); //0: don't split
}

template <typename DType> void DEventWriterROOT::Create_Branch_FundamentalArray(TTree* locTree, string locBranchName, string locArraySizeString, unsigned int locInitialSize) const
{
	string locTypeString = DROOTTypeString<DType>::GetTypeString();
	string locArrayName = locBranchName + string("[") + locArraySizeString + string("]/") + locTypeString;
	locTree->Branch(locBranchName.c_str(), new DType[locInitialSize], locArrayName.c_str());
	Get_FundamentalArraySizeMap()[locTree][locBranchName] = locInitialSize;
}

inline void DEventWriterROOT::Create_Branch_ClonesArray(TTree* locTree, string locBranchName, string locClassName, unsigned int locSize) const
{
	Get_ClonesArrayMap()[locTree].insert(pair<string, TClonesArray*>(locBranchName, new TClonesArray(locClassName.c_str(), locSize)));
	locTree->Branch(locBranchName.c_str(), &(Get_ClonesArrayMap()[locTree][locBranchName]), 32000, 0); //0: don't split
}

//BRANCH FILLING: //with separate particle & variable names (from which the branch name is made)
template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Fill_FundamentalData<DType>(locTree, locBranchName, locValue);
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex) const
{
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Fill_FundamentalData<DType>(locTree, locBranchName, locValue, locArrayIndex);
}

template <typename DType> void DEventWriterROOT::Fill_ClonesData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject, unsigned int locArrayIndex) const
{
	//only call for objects inheriting from TObject*!!!
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Fill_ClonesData<DType>(locTree, locBranchName, locObject, locArrayIndex);
}

template <typename DType> void DEventWriterROOT::Fill_TObjectData(TTree* locTree, string locParticleBranchName, string locVariableName, DType& locObject) const
{
	//only call for objects inheriting from TObject*!!!
	string locBranchName = Build_BranchName(locParticleBranchName, locVariableName);
	Fill_TObjectData<DType>(locTree, locBranchName, locObject);
}

//BRANCH FILLING: //with the full branch name
template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue) const
{
	DType* locBranchPointer = (DType*)locTree->GetBranch(locBranchName.c_str())->GetAddress();
	*locBranchPointer = locValue;
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locBranchName, DType locValue, unsigned int locArrayIndex) const
{
	//create a new, larger array if the current one is too small
	unsigned int& locCurrentArraySize = Get_FundamentalArraySizeMap()[locTree][locBranchName];
	DType* locBranchPointer = (DType*)locTree->GetBranch(locBranchName.c_str())->GetAddress();
	if(locArrayIndex >= locCurrentArraySize)
	{
		DType* locOldBranchAddress = locBranchPointer;
		locTree->SetBranchAddress(locBranchName.c_str(), new DType[locArrayIndex + 1]);
		locBranchPointer = (DType*)locTree->GetBranch(locBranchName.c_str())->GetAddress();
		//copy the old contents into the new array
		for(unsigned int loc_i = 0; loc_i < locCurrentArraySize; ++loc_i)
			locBranchPointer[loc_i] = locOldBranchAddress[loc_i];
		delete[] locOldBranchAddress;
		locCurrentArraySize = locArrayIndex + 1;
	}

	locBranchPointer[locArrayIndex] = locValue;
}

template <typename DType> void DEventWriterROOT::Fill_ClonesData(TTree* locTree, string locBranchName, DType& locObject, unsigned int locArrayIndex) const
{
	//only call for objects inheriting from TObject*!!!
	TClonesArray* locClonesArray = Get_ClonesArrayMap()[locTree][locBranchName];
	DType* locConstructedObject = (DType*)locClonesArray->ConstructedAt(locArrayIndex);
	*locConstructedObject = locObject;
}

template <typename DType> void DEventWriterROOT::Fill_TObjectData(TTree* locTree, string locBranchName, DType& locObject) const
{
	//only call for objects inheriting from TObject*!!!
	DType* locDType = (DType*)Get_TObjectMap()[locTree][locBranchName];
	*locDType = locObject;
}

inline string DEventWriterROOT::Convert_ToBranchName(string locInputName) const
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

inline string DEventWriterROOT::Build_BranchName(string locParticleBranchName, string locVariableName) const
{
	return ((locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName);
}

#endif //_DEventWriterROOT_
