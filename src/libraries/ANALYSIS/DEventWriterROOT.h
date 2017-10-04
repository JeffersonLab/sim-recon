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

#include "JANA/JApplication.h"
#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"

#include "TRIGGER/DTrigger.h"
#include "BCAL/DBCALShower.h"
#include "FCAL/DFCALShower.h"
#include "TRACKING/DMCThrown.h"
#include "TRACKING/DTrackTimeBased.h"

#include "PID/DVertex.h"
#include "PID/DChargedTrack.h"
#include "PID/DBeamPhoton.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralShower.h"
#include "PID/DEventRFBunch.h"
#include "PID/DMCReaction.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisResults.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DMCThrownMatching.h"
#include "ANALYSIS/DCutActions.h"
#include "ANALYSIS/DTreeInterface.h"
#include "ANALYSIS/DReactionVertexInfo.h"

using namespace std;
using namespace jana;

class DEventWriterROOT : public JObject
{
	public:
		JOBJECT_PUBLIC(DEventWriterROOT);

		virtual ~DEventWriterROOT(void);
		void Initialize(JEventLoop* locEventLoop);

		void Create_ThrownTree(JEventLoop* locEventLoop, string locOutputFileName) const;

		void Fill_DataTrees(JEventLoop* locEventLoop, string locDReactionTag) const; //fills all from this factory tag
		void Fill_DataTree(JEventLoop* locEventLoop, const DReaction* locReaction, deque<const DParticleCombo*>& locParticleCombos) const;
		void Fill_ThrownTree(JEventLoop* locEventLoop) const;

	protected:

		//CUSTOM FUNCTIONS: //Inherit from this class and write custom code in these functions
			//DO NOT: Write any code that requires a lock of ANY KIND. No reading calibration constants, accessing gParams, etc. This can cause deadlock.
		virtual void Create_CustomBranches_ThrownTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop) const{};
		virtual void Fill_CustomBranches_ThrownTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns) const{};
		virtual void Create_CustomBranches_DataTree(DTreeBranchRegister& locBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction, bool locIsMCDataFlag) const{};
		virtual void Fill_CustomBranches_DataTree(DTreeFillData* locTreeFillData, JEventLoop* locEventLoop, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const DMCThrownMatching* locMCThrownMatching, const DDetectorMatches* locDetectorMatches,
				const vector<const DBeamPhoton*>& locBeamPhotons, const vector<const DChargedTrackHypothesis*>& locChargedHypos,
				const vector<const DNeutralParticleHypothesis*>& locNeutralHypos, const deque<const DParticleCombo*>& locParticleCombos) const{};

		//UTILITY FUNCTIONS
		string Convert_ToBranchName(string locInputName) const;
		string Build_BranchName(string locParticleBranchName, string locVariableName) const;
		ULong64_t Calc_ParticleMultiplexID(Particle_t locPID) const;
		void Get_DecayProductNames(const DReaction* locReaction, size_t locReactionStepIndex, TMap* locPositionToNameMap, TList*& locDecayProductNames, deque<size_t>& locSavedSteps) const;

		const DAnalysisUtilities* dAnalysisUtilities;

	private:

		unsigned int dInitNumThrownArraySize;
		unsigned int dInitNumBeamArraySize;
		unsigned int dInitNumTrackArraySize;
		unsigned int dInitNumNeutralArraySize;
		unsigned int dInitNumComboArraySize;

		double dTargetCenterZ;

		//DEFAULT ACTIONS LISTED SEPARATELY FROM CUSTOM (in case in derived class user does something bizarre)
		map<const DReaction*, DCutAction_ThrownTopology*> dCutActionMap_ThrownTopology;
		map<const DReaction*, DCutAction_TrueCombo*> dCutActionMap_TrueCombo;
		map<const DReaction*, DCutAction_BDTSignalCombo*> dCutActionMap_BDTSignalCombo;

		//add in future: let user execute custom actions (outside of lock): user adds and initializes actions in derived-writer constructor
		//map<const DReaction*, map<string, map<const DParticleCombo*, bool> > > dCustomActionResultsMap; //string is action name

		/****************************************************************************************************************************************/

		//TREE INTERFACES, FILL OBJECTS
		//The non-thrown objects are created during the constructor, and thus the maps can remain const
		//The thrown objects are created later by the user (so they can specify file name), when the object is const, so they are declared mutable
		mutable DTreeInterface* dThrownTreeInterface;
		mutable DTreeFillData dThrownTreeFillData;
		map<const DReaction*, DTreeInterface*> dTreeInterfaceMap;
		map<const DReaction*, DTreeFillData*> dTreeFillDataMap;

		map<const DReaction*, const DReactionVertexInfo*> dVertexInfoMap;

		//TREE CREATION:
		void Create_DataTree(const DReaction* locReaction, JEventLoop* locEventLoop, bool locIsMCDataFlag);
		TMap* Create_UserInfoMaps(DTreeBranchRegister& locTreeBranchRegister, JEventLoop* locEventLoop, const DReaction* locReaction) const;
		void Create_UserTargetInfo(DTreeBranchRegister& locTreeBranchRegister, Particle_t locTargetPID) const;
		void Create_Branches_Thrown(DTreeBranchRegister& locTreeBranchRegister, bool locIsOnlyThrownFlag) const;

		//TREE CREATION: PARTICLE INFO
		void Create_Branches_ThrownParticles(DTreeBranchRegister& locTreeBranchRegister, bool locIsOnlyThrownFlag) const;
		void Create_Branches_Beam(DTreeBranchRegister& locTreeBranchRegister, bool locIsMCDataFlag) const;
		void Create_Branches_NeutralHypotheses(DTreeBranchRegister& locTreeBranchRegister, bool locIsMCDataFlag) const;
		void Create_Branches_ChargedHypotheses(DTreeBranchRegister& locTreeBranchRegister, bool locIsMCDataFlag) const;

		//TREE CREATION: COMBO INFO
			//TMap is locPositionToNameMap
		void Create_Branches_Combo(DTreeBranchRegister& locTreeBranchRegister, const DReaction* locReaction, bool locIsMCDataFlag, TMap* locPositionToNameMap) const;
		void Create_Branches_BeamComboParticle(DTreeBranchRegister& locTreeBranchRegister, Particle_t locBeamPID, DKinFitType locKinFitType) const;
		void Create_Branches_ComboTrack(DTreeBranchRegister& locTreeBranchRegister, string locParticleBranchName, DKinFitType locKinFitType) const;
		void Create_Branches_ComboNeutral(DTreeBranchRegister& locTreeBranchRegister, string locParticleBranchName, DKinFitType locKinFitType) const;

		//TREE FILLING: THROWN INFO
		void Compute_ThrownPIDInfo(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying,
				ULong64_t& locNumPIDThrown_FinalState, ULong64_t& locPIDThrown_Decaying) const;
		void Group_ThrownParticles(const vector<const DMCThrown*>& locMCThrowns_FinalState, const vector<const DMCThrown*>& locMCThrowns_Decaying,
				vector<const DMCThrown*>& locMCThrownsToSave, map<const DMCThrown*, unsigned int>& locThrownIndexMap) const;
		void Fill_ThrownInfo(DTreeFillData* locTreeFillData, const DMCReaction* locMCReaction, const vector<const DMCThrown*>& locMCThrowns,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, ULong64_t locNumPIDThrown_FinalState, ULong64_t locPIDThrown_Decaying,
				const DMCThrownMatching* locMCThrownMatching = NULL) const;
		void Fill_ThrownParticleData(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DMCThrown* locMCThrown, 
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DMCThrownMatching* locMCThrownMatching) const;

		//TREE FILLING: GET HYPOTHESES/BEAM
		vector<const DBeamPhoton*> Get_BeamPhotons(const deque<const DParticleCombo*>& locParticleCombos) const;
		vector<const DChargedTrackHypothesis*> Get_ChargedHypotheses(JEventLoop* locEventLoop) const;
		vector<const DChargedTrackHypothesis*> Get_ChargedHypotheses_Used(JEventLoop* locEventLoop, const DReaction* locReaction, const deque<const DParticleCombo*>& locParticleCombos) const;
		vector<const DNeutralParticleHypothesis*> Get_NeutralHypotheses(JEventLoop* locEventLoop, const set<Particle_t>& locReactionPIDs) const;
		vector<const DNeutralParticleHypothesis*> Get_NeutralHypotheses_Used(JEventLoop* locEventLoop, const DReaction* locReaction, const set<Particle_t>& locReactionPIDs, const deque<const DParticleCombo*>& locParticleCombos) const;

		//TREE FILLING: INDEPENDENT PARTICLES
		void Fill_BeamData(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DBeamPhoton* locBeamPhoton, const DVertex* locVertex, const DMCThrownMatching* locMCThrownMatching) const;
		void Fill_ChargedHypo(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DChargedTrackHypothesis* locChargedTrackHypothesis, const DMCThrownMatching* locMCThrownMatching,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const;
		void Fill_NeutralHypo(DTreeFillData* locTreeFillData, unsigned int locArrayIndex, const DNeutralParticleHypothesis* locPhotonHypothesis, const DMCThrownMatching* locMCThrownMatching,
				const map<const DMCThrown*, unsigned int>& locThrownIndexMap, const DDetectorMatches* locDetectorMatches) const;

		//TREE FILLING: COMBO
		void Fill_ComboData(DTreeFillData* locTreeFillData, const DReaction* locReaction, const DParticleCombo* locParticleCombo, unsigned int locComboIndex, const map<pair<oid_t, Particle_t>, size_t>& locObjectToArrayIndexMap) const;
		void Fill_ComboStepData(DTreeFillData* locTreeFillData, const DReaction* locReaction, const DParticleCombo* locParticleCombo, unsigned int locStepIndex, unsigned int locComboIndex,
				DKinFitType locKinFitType, const map<pair<oid_t, Particle_t>, size_t>& locObjectToArrayIndexMap) const;

		//TREE FILLING: COMBO PARTICLES
		void Fill_ComboBeamData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, const DBeamPhoton* locBeamPhoton, size_t locBeamIndex, DKinFitType locKinFitType) const;
		void Fill_ComboChargedData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, string locParticleBranchName, const DChargedTrackHypothesis* locMeasuredChargedHypo,
				const DChargedTrackHypothesis* locChargedHypo, size_t locChargedIndex, DKinFitType locKinFitType) const;
		void Fill_ComboNeutralData(DTreeFillData* locTreeFillData, unsigned int locComboIndex, string locParticleBranchName, const DNeutralParticleHypothesis* locMeasuredNeutralHypo,
				const DNeutralParticleHypothesis* locNeutralHypo, size_t locNeutralIndex, DKinFitType locKinFitType) const;
};

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
