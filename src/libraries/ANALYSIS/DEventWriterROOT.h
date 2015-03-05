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
			//DO NOT: Acquire/release the ROOT lock.  It is already acquired prior to entry into these functions
			//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
			//DO NOT: Call TTree::Fill().  This will be called after calling the custom fill functions.
			//DO: Use the inherited functions for creating/filling branches.  They will make your life MUCH easier.
		virtual void Create_CustomBranches_DataTree(TTree* locTree, const DReaction* locReaction, bool locIsMCDataFlag) const{};
		virtual void Create_CustomBranches_ThrownTree(TTree* locTree) const{};
		virtual void Fill_CustomBranches_DataTree(TTree* locTree, JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo) const{};
		virtual void Fill_CustomBranches_ThrownTree(TTree* locTree, JEventLoop* locEventLoop) const{};

		//UTILITY FUNCTIONS
		string Convert_ToBranchName(string locInputName) const;
		ULong64_t Calc_ParticleMultiplexID(Particle_t locPID) const;
		void Get_DecayProductNames(const DReaction* locReaction, size_t locReactionStepIndex, TMap* locPositionToNameMap, TList*& locDecayProductNames, deque<size_t>& locSavedSteps) const;

		//BRANCH CREATION:
		template <typename DType> string Create_Branch_Fundamental(TTree* locTree, string locParticleBranchName, string locVariableName, string locTypeString) const;
		template <typename DType> string Create_Branch_NoSplitTObject(TTree* locTree, string locParticleBranchName, string locVariableName, map<string, TObject*>& locTObjectMap) const;
		template <typename DType> string Create_Branch_FundamentalArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locArraySizeString, unsigned int locMinimumSize, string locTypeString) const;
		string Create_Branch_ClonesArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locClassName, unsigned int locSize, map<string, TClonesArray*>& locClonesArrayMap) const;

		//TREE FILLING:
		template <typename DType> DType* Get_BranchAddress(TTree* locTree, string locBranchName, unsigned int locMinimumSize, unsigned int locCurrentSize) const;
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locVariableName, DType locValue) const;
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue) const;
		template <typename DType> void Fill_ClonesData(TTree* locTree, string locParticleBranchName, string locVariableName, DType* locObject, unsigned int locArrayIndex, const map<string, TClonesArray*>& locClonesArrayMap) const;
		template <typename DType> void Fill_TObjectData(TTree* locTree, string locParticleBranchName, string locVariableName, DType* locObject, const map<string, TObject*>& locTObjectMap) const;

		//calls the next-below method with locMinArraySize = locArrayIndex
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex, unsigned int locCurrentArraySize) const;
		//locMinArraySize: the minimum array size needed for this event: if it is less than the current size, a new, larger array will be created
			//this method is useful for auto-expanding the array to the size you know you need, rather than calling new for each new array index
		template <typename DType> void Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex, unsigned int locMinArraySize, unsigned int locCurrentArraySize) const;

		const DAnalysisUtilities* dAnalysisUtilities;

	private:
		DEventWriterROOT(void){}; //don't allow default constructor

		void Get_Reactions(jana::JEventLoop* locEventLoop, vector<const DReaction*>& locReactions) const;

		//TREE CREATION:
		void Create_DataTree(const DReaction* locReaction, bool locIsMCDataFlag) const;
		void Create_Branches_FinalStateParticle(TTree* locTree, string locParticleBranchName, bool locIsChargedFlag, bool locKinFitFlag, bool locIsMCDataFlag) const;
		void Create_Branches_Beam(TTree* locTree, string locParticleBranchName, bool locKinFitFlag) const;
		void Create_Branches_UnusedParticle(TTree* locTree, string locParticleBranchName, string locArraySizeString, bool locIsMCDataFlag) const;
		void Create_Branches_ThrownParticle(TTree* locTree, string locParticleBranchName, string locArraySizeString, bool locIsOnlyThrownFlag) const;

		//TREE FILLING:
		void Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DMCThrown* locMCThrown, map<const DMCThrown*, unsigned int> locThrownObjectIDMap) const;
		void Fill_ThrownParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DMCThrown* locMCThrown, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, const map<const DNeutralShower*, int>& locShowerToIDMap) const;
		void Fill_UnusedParticleData(TTree* locTree, unsigned int locArrayIndex, unsigned int locMinArraySize, const DKinematicData* locKinematicData, const DEventRFBunch* locEventRFBunch, const map<const DNeutralShower*, int>& locShowerToIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DDetectorMatches* locDetectorMatches) const;
		void Fill_BeamParticleData(TTree* locTree, string locParticleBranchName, const DKinematicData* locKinematicData, const DKinematicData* locKinematicData_Measured, const map<const DBeamPhoton*, int>& locBeamToIDMap) const;
		void Fill_ParticleData(bool locKinFitFlag, TTree* locTree, string locParticleBranchName, const DKinematicData* locKinematicData, const DKinematicData* locKinematicData_Measured, const DEventRFBunch* locEventRFBunch, const map<const DNeutralShower*, int>& locShowerToIDMap, const DMCThrownMatching* locMCThrownMatching, double locMinThrownMatchFOM, map<const DMCThrown*, unsigned int> locThrownObjectIDMap, const DDetectorMatches* locDetectorMatches) const;
};

template <typename DType> string DEventWriterROOT::Create_Branch_Fundamental(TTree* locTree, string locParticleBranchName, string locVariableName, string locTypeString) const
{
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	string locTypeName = locBranchName + string("/") + locTypeString;
	locTree->Branch(locBranchName.c_str(), new DType(), locTypeName.c_str());
	return locBranchName;
}

template <typename DType> string DEventWriterROOT::Create_Branch_NoSplitTObject(TTree* locTree, string locParticleBranchName, string locVariableName, map<string, TObject*>& locTObjectMap) const
{
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	locTObjectMap.insert(pair<string, TObject*>(locBranchName, (TObject*)(new DType())));
	locTree->Branch<DType>(locBranchName.c_str(), (DType**)&(locTObjectMap[locBranchName]), 32000, 0); //0: don't split
	return locBranchName;
}

template <typename DType> string DEventWriterROOT::Create_Branch_FundamentalArray(TTree* locTree, string locParticleBranchName, string locVariableName, string locArraySizeString, unsigned int locMinimumSize, string locTypeString) const
{
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	string locArrayName = locBranchName + string("[") + locArraySizeString + string("]/") + locTypeString;
	locTree->Branch(locBranchName.c_str(), new DType[locMinimumSize], locArrayName.c_str());
	return locBranchName;
}

template <typename DType> DType* DEventWriterROOT::Get_BranchAddress(TTree* locTree, string locBranchName, unsigned int locMinimumSize, unsigned int locCurrentSize) const
{
	//only works for fundamental types!!!
	if(locMinimumSize > locCurrentSize) //creates a new array if it is too small
	{
		delete[] (DType*)locTree->GetBranch(locBranchName.c_str())->GetAddress();
		locTree->SetBranchAddress(locBranchName.c_str(), new DType[locMinimumSize]);
	}
	return (DType*)locTree->GetBranch(locBranchName.c_str())->GetAddress();
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue) const
{
	//only call if the variable is NOT an array!!
	Fill_FundamentalData<DType>(locTree, locParticleBranchName, locVariableName, locValue, 0, 1, 1);
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locVariableName, DType locValue) const
{
	//only call if the variable is NOT an array!!
	Fill_FundamentalData<DType>(locTree, "", locVariableName, locValue, 0, 1, 1);
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex, unsigned int locCurrentArraySize) const
{
	Fill_FundamentalData(locTree, locParticleBranchName, locVariableName, locValue, locArrayIndex, locArrayIndex, locCurrentArraySize);
}

template <typename DType> void DEventWriterROOT::Fill_FundamentalData(TTree* locTree, string locParticleBranchName, string locVariableName, DType locValue, unsigned int locArrayIndex, unsigned int locMinArraySize, unsigned int locCurrentArraySize) const
{
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	if(locArrayIndex >= locMinArraySize)
		locMinArraySize = locArrayIndex + 1; //in case a user screws this up
	DType* locBranchPointer = Get_BranchAddress<DType>(locTree, locBranchName, locMinArraySize, locCurrentArraySize);
	if(locMinArraySize == 1)
		*locBranchPointer = locValue;
	else
		locBranchPointer[locArrayIndex] = locValue;
}

template <typename DType> void DEventWriterROOT::Fill_ClonesData(TTree* locTree, string locParticleBranchName, string locVariableName, DType* locObject, unsigned int locArrayIndex, const map<string, TClonesArray*>& locClonesArrayMap) const
{
	//only call for objects inheriting from TObject*!!!
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	TClonesArray* locClonesArray = locClonesArrayMap.find(locBranchName)->second;
	DType* locConstructedObject = (DType*)locClonesArray->ConstructedAt(locArrayIndex);
	*locConstructedObject = *locObject;
}

template <typename DType> void DEventWriterROOT::Fill_TObjectData(TTree* locTree, string locParticleBranchName, string locVariableName, DType* locObject, const map<string, TObject*>& locTObjectMap) const
{
	//only call for objects inheriting from TObject*!!!
	string locBranchName = (locParticleBranchName != "") ? locParticleBranchName + string("__") + locVariableName : locVariableName;
	DType* locDType = (DType*)locTObjectMap.find(locBranchName)->second;
	*locDType = *locObject;
}

#endif //_DEventWriterROOT_
