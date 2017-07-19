#ifndef DSourceComboer_h
#define DSourceComboer_h

#include <map>
#include <set>
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <stack>

#include "TF1.h"

#include "JANA/JObject.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "SplitString.h"
#include "DANA/DApplication.h"
#include "HDGEOMETRY/DGeometry.h"
#include "EVENTSTORE/DESSkimData.h"

#include "PID/DNeutralShower.h"
#include "PID/DKinematicData.h"
#include "PID/DEventRFBunch.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DSourceCombo.h"
#include "ANALYSIS/DReactionStepVertexInfo.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "DResourcePool.h"

#include "ANALYSIS/DSourceComboVertexer.h"
#include "ANALYSIS/DSourceComboP4Handler.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DParticleComboCreator.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

/****************************************************** DEFINE LAMBDAS, USING STATEMENTS *******************************************************/

struct DCompare_SourceComboInfos{
	bool operator()(const DSourceComboInfo* lhs, const DSourceComboInfo* rhs) const{return *lhs < *rhs;}
};

/********************************************************** DEFINE USING STATEMENTS ***********************************************************/

//DEFINE USING STATEMENTS
using DCombosByBeamBunch = map<vector<int>, vector<const DSourceCombo*>>;
using DSourceCombosByBeamBunchByUse = map<DSourceComboUse, DCombosByBeamBunch>;
//The DSourceCombosByUse_Large type uses a vector to pointer so that the combos can be easily copied and reused for another use
//e.g. when you can't place a mass cut yet: 2 different uses, identical combos: far faster to just copy the pointer to the large vector
using DSourceCombosByUse_Large = map<DSourceComboUse, vector<const DSourceCombo*>*>;
using DCombosByReaction = unordered_map<const DReaction*, vector<const DParticleCombo*>>;

/************************************************************** DEFINE CLASSES ***************************************************************/

class DSourceComboer : public JObject
{
	enum ComboingStage_t
	{
		d_ChargedStage,
		d_MixedStage_ZIndependent,
		d_MixedStage
	};

	enum class DConstructionStage
	{
		Min_Particles = 0,
		Max_Particles,
		In_Skim,
		Charged_Combos,
		Charged_RFBunch,
		Full_Combos,
		Neutral_RFBunch,
		NoVertex_RFBunch,
		HeavyNeutral_IM,
		Beam_Combos,
		MMVertex_Timing,
		MMVertex_IMCuts,
		Reaction_BeamRFCuts,
		Missing_Mass
	};

	using DConstructionStageType = std::underlying_type<DConstructionStage>::type;

	public:

		DSourceComboer(void) = delete;
		DSourceComboer(JEventLoop* locEventLoop);
		~DSourceComboer(void);

		//RESET
		void Reset_NewEvent(JEventLoop* locEventLoop);

		//BUILD COMBOS (what should be called from the outside to do all of the work)
		DCombosByReaction Build_ParticleCombos(const DReactionVertexInfo* locReactionVertexInfo);

		//Get combo characteristics
		Charge_t Get_ChargeContent(const DSourceComboInfo* locSourceComboInfo) const{return dComboInfoChargeContent.find(locSourceComboInfo)->second;}
		bool Get_HasMassiveNeutrals(const DSourceComboInfo* locSourceComboInfo) const{return (dComboInfosWithMassiveNeutrals.find(locSourceComboInfo) != dComboInfosWithMassiveNeutrals.end());}

		//Combo utility functions
		const DSourceCombo* Get_StepSourceCombo(const DReaction* locReaction, size_t locDesiredStepIndex, const DSourceCombo* locVertexPrimaryCombo, size_t locVertexPrimaryStepIndex = 0) const;
		const DSourceCombo* Get_VertexPrimaryCombo(const DSourceCombo* locReactionCombo, const DReactionStepVertexInfo* locStepVertexInfo) const;
		const DSourceCombo* Get_VertexPrimaryCombo(const DSourceCombo* locReactionCombo, const DReactionStepVertexInfo* locStepVertexInfo);

		//Get combo uses
		DSourceComboUse Get_SourceComboUse(const DReactionStepVertexInfo* locStepVertexInfo) const{return dSourceComboUseReactionMap.find(locStepVertexInfo)->second;};
		DSourceComboUse Get_SourceComboUse(const DReaction* locReaction, size_t locStepIndex) const{return dSourceComboUseReactionStepMap.find(locReaction)->second.find(locStepIndex)->second;};
		DSourceComboUse Get_PrimaryComboUse(const DReactionVertexInfo* locReactionVertexInfo) const{return Get_SourceComboUse(locReactionVertexInfo->Get_StepVertexInfo(0));};

		DParticleComboCreator* Get_ParticleComboCreator(void) const{return dParticleComboCreator;}
		void Print_NumCombosByUse(void);

	private:

		/********************************************************** DECLARE MEMBER FUNCTIONS ***********************************************************/

		//SETUP
		void Setup_NeutralShowers(JEventLoop* locEventLoop);

		//INITIAL CHECKS
		bool Check_Reactions(vector<const DReaction*>& locReactions);
		bool Check_NumParticles(const DReaction* locReaction);
		bool Check_Skims(const DReaction* locReaction) const;

		//PARTICLE CUTS
		bool Cut_dEdxAndEOverP(const DChargedTrackHypothesis* locHypo);
		bool Cut_dEdx(Particle_t locPID, DetectorSystem_t locSystem, double locP, double locdEdx);
		bool Cut_EOverP(Particle_t locPID, DetectorSystem_t locSystem, double locP, double locEOverP);
		void Fill_CutHistograms(void);
		void Fill_SurvivalHistograms(void);

		//CREATE PHOTON COMBO INFOS & USES
		void Create_SourceComboInfos(const DReactionVertexInfo* locReactionVertexInfo);
		DSourceComboUse Create_ZDependentSourceComboUses(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionChargedCombo);
		DSourceComboUse Build_NewZDependentUse(const DReaction* locReaction, size_t locStepIndex, signed char locVertexZBin, const DSourceComboUse& locOrigUse, const unordered_map<size_t, DSourceComboUse>& locCreatedUseMap);

		//CREATE PHOTON COMBO INFOS & USES: UTILITY METHODS
		map<Particle_t, unsigned char> Build_ParticleMap(const DReaction* locReaction, size_t locStepIndex, Charge_t locCharge) const;
		pair<bool, map<DSourceComboUse, unsigned char>> Get_FinalStateDecayingComboUses(const DReaction* locReaction, size_t locStepIndex, const map<size_t, DSourceComboUse>& locStepComboUseMap) const;
		DSourceComboUse Make_ComboUse(Particle_t locInitPID, const map<Particle_t, unsigned char>& locNumParticles, const map<DSourceComboUse, unsigned char>& locFurtherDecays);
		const DSourceComboInfo* MakeOrGet_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays, unsigned char locNumTabs);
		const DSourceComboInfo* GetOrMake_SourceComboInfo(const vector<pair<Particle_t, unsigned char>>& locNumParticles, const vector<pair<DSourceComboUse, unsigned char>>& locFurtherDecays, unsigned char locNumTabs);

		//CREATE COMBOS
		void Combo_WithNeutralsAndBeam(const vector<const DReaction*>& locReactions, const DReactionVertexInfo* locReactionVertexInfo, const DSourceComboUse& locPrimaryComboUse, const DSourceCombo* locReactionChargedCombo, const vector<int>& locBeamBunches_Charged, DCombosByReaction& locOutputComboMap);
		void Combo_WithBeam(const vector<const DReaction*>& locReactions, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locReactionFullCombo, int locRFBunch, DCombosByReaction& locOutputComboMap);
		const DParticleCombo* Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle);

		//CREATE SOURCE COMBOS - GENERAL METHODS
		void Create_SourceCombos(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);
		void Create_SourceCombos_Unknown(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);

		//COMBO VERTICALLY METHODS
		//Note that vertical comboing always takes place at the same vertex-z
		void Combo_Vertically_AllDecays(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);
		void Combo_Vertically_NDecays(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, const DSourceComboUse& locSourceComboDecayUse, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);
		void Combo_Vertically_AllParticles(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, unsigned char locNumTabs);
		void Combo_Vertically_NParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locNMinus1ComboUse, ComboingStage_t locComboingStage, unsigned char locNumTabs);

		//COMBO HORIZONTALLY ORGANIZATION METHODS
		void Combo_Horizontally_All(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);
		void Combo_Horizontally_AddDecay(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locComboUseAllBut1, const DSourceComboUse& locComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs);
		void Combo_Horizontally_AddParticles(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locComboUseAllBut1, const pair<Particle_t, unsigned char>& locParticlePairToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs);

		//COMBO HORIZONTALLY COMBOING METHODS
		void Create_Combo_OneParticle(const DSourceComboUse& locComboUseToCreate, ComboingStage_t locComboingStage, unsigned char locNumTabs);
		void Combo_Horizontally_AddCombo(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, const DSourceComboUse& locSourceComboUseToAdd, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, bool locExpandAllBut1Flag, unsigned char locNumTabs);
		void Combo_Horizontally_AddParticle(const DSourceComboUse& locComboUseToCreate, const DSourceComboUse& locAllBut1ComboUse, Particle_t locPID, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_Presiding, unsigned char locNumTabs);

		//BUILD/RETRIEVE RESUME-AT ITERATORS
		void Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles);
		void Build_ComboIndices(const DSourceComboUse& locSourceComboUse, const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage);
		vector<const JObject*>::const_iterator Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const;
		size_t Get_ResumeAtIndex_Combos(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locPreviousCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage) const;

		//GET POTENTIAL PARTICLES & COMBOS FOR COMBOING
		const vector<const DSourceCombo*>& Get_CombosForComboing(const DSourceComboUse& locComboUse, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches, const DSourceCombo* locChargedCombo_WithPrevious);
		const vector<const DSourceCombo*>& Get_CombosByBeamBunch(const DSourceComboUse& locComboUse, DCombosByBeamBunch& locCombosByBunch, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage);

		//REGISTER VALID RF BUNCHES
		void Register_ValidRFBunches(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locSourceCombo, const vector<int>& locRFBunches, ComboingStage_t locComboingStage, const DSourceCombo* locChargedCombo_WithNow);

		//PARTICLE UTILITY FUNCTIONS
		const vector<const JObject*>& Get_ParticlesForComboing(Particle_t locPID, ComboingStage_t locComboingStage, const vector<int>& locBeamBunches = {}, signed char locVertexZBin = 0);
		const vector<const JObject*>& Get_ShowersByBeamBunch(const vector<int>& locBeamBunches, DPhotonShowersByBeamBunch& locShowersByBunch);
		shared_ptr<const DKinematicData> Create_KinematicData(const DNeutralShower* locNeutralShower, const DVector3& locVertex) const;
		bool Get_IsComboingZIndependent(const JObject* locObject, Particle_t locPID) const;

		//COMBO UTILITY FUNCTIONS
		DSourceCombosByUse_Large& Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		DSourceCombosByBeamBunchByUse& Get_SourceCombosByBeamBunchByUse(Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo = nullptr);
		void Copy_ZIndependentMixedResults(const DSourceComboUse& locComboUseToCreate, const DSourceCombo* locChargedCombo_WithNow);
		const DSourceCombo* Get_ChargedCombo_WithNow(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboInfo* locToCreateComboInfo, ComboingStage_t locComboingStage) const;
		const DSourceCombo* Get_NextChargedCombo(const DSourceCombo* locChargedCombo_Presiding, const DSourceComboUse& locNextComboUse, ComboingStage_t locComboingStage, bool locGetPresidingFlag, size_t locInstance) const;
		bool Get_PromoteFlag(Particle_t locDecayPID_UseToCheck, const DSourceComboInfo* locComboInfo_UseToCreate, const DSourceComboInfo* locComboInfo_UseToCheck) const;
		const DSourceCombo* Find_Combo_AtThisStep(const DSourceCombo* locSourceCombo, DSourceComboUse locUseToFind, size_t locDecayInstanceIndex) const;

		//GET RESOURCES
		DSourceCombo* Get_SourceComboResource(void);
		vector<const DSourceCombo*>* Get_SourceComboVectorResource(void);

		/************************************************************** DEFINE MEMBERS ***************************************************************/

		//CONTROL INFORMATION
		uint64_t dEventNumber = 0; //re-setup on new events
		string dShowerSelectionTag = "PreSelect";
		size_t dDebugLevel = 0;

		//EXPERIMENT INFORMATION
		DVector3 dTargetCenter;

		//COMMAND LINE CUTS
		pair<bool, size_t> dNumPlusMinusRFBunches = std::make_pair(false, 0); //by default use DReaction cut //only use this if set on command line
		unordered_map<const DReaction*, size_t> dRFBunchCutsByReaction;
		unordered_map<const DReactionVertexInfo*, size_t> dMaxRFBunchCuts;

		//HANDLERS AND VERTEXERS
		DSourceComboVertexer* dSourceComboVertexer;
		DSourceComboP4Handler* dSourceComboP4Handler;
		DSourceComboTimeHandler* dSourceComboTimeHandler;
		DParticleComboCreator* dParticleComboCreator;

		//SOURCE COMBO INFOS: CREATED ONCE DURING DSOURCECOMBOER OBJECT CONSTRUCTION
			//with some exceptions (specific vertex-z, etc.)
		//want to make sure we only have one of each type: suggests using a set
		//however, after the first few events, almost all of these have already been created: vector has faster lookup time
		//therefore, use the set when creating the objects during construction, but then move the results into the vector and keep it sorted
		set<const DSourceComboInfo*, DCompare_SourceComboInfos> dSourceComboInfoSet;
		vector<const DSourceComboInfo*> dSourceComboInfos;
		unordered_map<const DSourceComboInfo*, Charge_t> dComboInfoChargeContent;
		unordered_set<const DSourceComboInfo*> dComboInfosWithMassiveNeutrals;
		//the rest
		unordered_map<const DReactionStepVertexInfo*, DSourceComboUse> dSourceComboUseReactionMap; //primary combo info (nullptr if none)
		//combo use -> step
		map<pair<const DReactionStepVertexInfo*, DSourceComboUse>, size_t> dSourceComboInfoStepMap; //size_t: step index
		//i need to go from step -> combo use
		unordered_map<const DReaction*, map<size_t, DSourceComboUse>> dSourceComboUseReactionStepMap; //primary combo info (nullptr if none)
		//with specific vertex-z's
		map<pair<const DReactionVertexInfo*, vector<signed char>>, DSourceComboUse> dSourceComboUseVertexZMap;
		map<DSourceComboUse, DSourceComboUse> dZDependentUseToIndependentMap; //from z-dependent -> z-independent

		//SKIM INFORMATION
		const DESSkimData* dESSkimData = nullptr;

		//PARTICLES
		size_t dMaxNumNeutrals = 20;
		map<Particle_t, vector<const JObject*>> dTracksByPID;
		size_t dNumChargedTracks;
		map<bool, vector<const JObject*>> dTracksByCharge; //true/false: positive/negative
		unordered_map<signed char, DPhotonShowersByBeamBunch> dShowersByBeamBunchByZBin; //char: zbin //for all showers: unknown z-bin, {} RF bunch

		//SOURCE COMBOS //vector: z-bin //if attempted and all failed, DSourceCombosByUse_Large vector will be empty
		size_t dInitialComboVectorCapacity = 100;
		DSourceCombosByUse_Large dSourceCombosByUse_Charged;
		unordered_map<const DSourceCombo*, DSourceCombosByUse_Large> dMixedCombosByUseByChargedCombo; //key: charged combo //value: contains mixed & neutral combos //neutral: key is nullptr
		//also, sort by which beam bunches they are valid for: that way when comboing, we can retrieve only the combos that can possibly match the input RF bunches
		unordered_map<const DSourceCombo*, DSourceCombosByBeamBunchByUse> dSourceCombosByBeamBunchByUse; //key: charged combo //value: contains mixed & neutral combos: key is nullptr
		map<pair<const DSourceCombo*, const DReactionStepVertexInfo*>, const DSourceCombo*> dVertexPrimaryComboMap; //first combo: reaction primary combo (can be charged or full!)

		//RESUME SEARCH ITERATORS
		//e.g. if a DSourceCombo is -> 2pi0, and we want to use it as a basis for building a combo of 3pi0s,
		//then this iterator points to the first pi0 in the DSourceCombosByUse_Large vector that we want to test
		//that way we save a lot of time, since we don't have to look for it again
		//they are useful when comboing VERTICALLY, but cannot be used when comboing HORIZONTALLY
			//e.g. when comboing a pi0 (with photons = A, D) with a single photon, the photon could be B, C, or E+: no single spot to resume at
		map<pair<const JObject*, vector<int>>, vector<const JObject*>::const_iterator> dResumeSearchAfterIterators_Particles; //vector<int>: RF bunches (empty for all)
		map<pair<const DSourceCombo*, DSourceComboUse>, map<vector<int>, size_t>> dResumeSearchAfterIndices_Combos; //char: zbin, size_t: index

		//RESUME-SEARCH-AFTER OBJECTS
		//The below resume-at vectors are only useful when comboing vertically, so PID-specific versions are not needed
		unordered_map<const DSourceCombo*, const JObject*> dResumeSearchAfterMap_Particles; //key: Combo containing N particles of the type in the value //value: the last of those N particles

		//VALID RF BUNCHES BY COMBO
		unordered_map<const DSourceCombo*, vector<int>> dValidRFBunches_ByCombo;

		//RESOURCE POOLS
		//Don't use these directly!  Use the Get_*Resource functions instead!!
		DResourcePool<DSourceCombo> dResourcePool_SourceCombo;
		DResourcePool<vector<const DSourceCombo*>> dResourcePool_SourceComboVector;
		//These are used to know what to recycle
		vector<DSourceCombo*> dCreatedCombos;
		vector<vector<const DSourceCombo*>*> dCreatedComboVectors;

		//Combo/event tracking
		map<const DReaction*, TH1*> dNumEventsSurvivedStageMap;
		map<const DReaction*, TH1*> dNumCombosSurvivedStageMap;
		map<const DReaction*, TH2*> dNumCombosSurvivedStage2DMap;
		map<const DReaction*, map<DConstructionStage, size_t>> dNumCombosSurvivedStageTracker; //index is for event stages!!!
		map<DSourceComboUse, size_t> dNumMixedCombosMap_Charged;
		map<DSourceComboUse, size_t> dNumMixedCombosMap_Mixed;

		//dE/dx
		map<Particle_t, map<DetectorSystem_t, pair<TF1*, TF1*>>> ddEdxCutMap; //pair: first is lower bound, second is upper bound
		map<Particle_t, map<DetectorSystem_t, vector<pair<double, double>>>> ddEdxValueMap; //pair: first is p, 2nd is dE/dx
		map<Particle_t, map<DetectorSystem_t, TH2*>> dHistMap_dEdx;

		//E/p
		map<Particle_t, map<DetectorSystem_t, TF1*>> dEOverPCutMap; //if lepton, select above function, else select below
		map<Particle_t, map<DetectorSystem_t, vector<pair<double, double>>>> dEOverPValueMap; //pair: first is p, 2nd is E/p
		map<Particle_t, map<DetectorSystem_t, TH2*>> dHistMap_EOverP;
};

/*********************************************************** INLINE MEMBER FUNCTION DEFINITIONS ************************************************************/


inline DSourceCombo* DSourceComboer::Get_SourceComboResource(void)
{
	auto locCombo = dResourcePool_SourceCombo.Get_Resource();
	locCombo->Reset();
	dCreatedCombos.push_back(locCombo);
	return locCombo;
}

inline vector<const DSourceCombo*>* DSourceComboer::Get_SourceComboVectorResource(void)
{
	auto locComboVector = dResourcePool_SourceComboVector.Get_Resource();
	locComboVector->clear();
	locComboVector->reserve(dInitialComboVectorCapacity);
	dCreatedComboVectors.push_back(locComboVector);
	return locComboVector;
}

inline bool DSourceComboer::Check_Skims(const DReaction* locReaction) const
{
	if(dESSkimData == nullptr)
		return true;

	string locReactionSkimString = locReaction->Get_EventStoreSkims();
	vector<string> locReactionSkimVector;
	SplitString(locReactionSkimString, locReactionSkimVector, ",");
	for(size_t loc_j = 0; loc_j < locReactionSkimVector.size(); ++loc_j)
	{
		if(!dESSkimData->Get_IsEventSkim(locReactionSkimVector[loc_j]))
			return false;
	}

	return true;
}

inline void DSourceComboer::Build_ParticleIterators(const vector<int>& locBeamBunches, const vector<const JObject*>& locParticles)
{
	for(vector<const JObject*>::const_iterator locIterator = locParticles.begin(); locIterator != locParticles.end(); ++locIterator)
		dResumeSearchAfterIterators_Particles.emplace(std::make_pair(*locIterator, locBeamBunches), locIterator);
}

inline void DSourceComboer::Build_ComboIndices(const DSourceComboUse& locSourceComboUse, const vector<int>& locBeamBunches, const vector<const DSourceCombo*>& locCombos, ComboingStage_t locComboingStage)
{
	for(size_t loc_i = 0; loc_i < locCombos.size(); ++loc_i)
	{
		if(dDebugLevel >= 20)
		{
			cout << "build resume iterators: vector address, combo, decay pid, zbin, index, bunches: " << &locCombos << ", " << locCombos[loc_i] << ", " << std::get<0>(locSourceComboUse) << ", " << int(std::get<1>(locSourceComboUse)) << ", " << loc_i << ", ";
			for(auto& locBunch : locBeamBunches)
				cout << locBunch << ", ";
			cout << endl;
		}
		dResumeSearchAfterIndices_Combos[std::make_pair(locCombos[loc_i], locSourceComboUse)].emplace(locBeamBunches, loc_i);
	}
}

inline vector<const JObject*>::const_iterator DSourceComboer::Get_ResumeAtIterator_Particles(const DSourceCombo* locSourceCombo, const vector<int>& locBeamBunches) const
{
	const auto& locPreviousObject = dResumeSearchAfterMap_Particles.find(locSourceCombo)->second;
	return std::next(dResumeSearchAfterIterators_Particles.find(std::make_pair(locPreviousObject, locBeamBunches))->second);
}

inline size_t DSourceComboer::Get_ResumeAtIndex_Combos(const DSourceComboUse& locSourceComboUse, const DSourceCombo* locPreviousCombo, const vector<int>& locBeamBunches, ComboingStage_t locComboingStage) const
{
	//Suppose you are comboing N pi0s together (let's say 3), and 6 different pi0 combos were reconstructed
	//So, you choose the first pi0, then want to loop over the remaining ones
	//To avoid duplication, you don't want to loop over ALL of the others, only the ones AFTER the first one
	//So, the input is the last pi0 that you've chosen for your combo
	//Then, just find the location where that pi0 is the possible-combos-vector (pi0 vector), and return that index + 1
	auto locSearchPair = std::make_pair(locPreviousCombo, locSourceComboUse);
	auto& locBunchIndexMap = dResumeSearchAfterIndices_Combos.find(locSearchPair)->second;
	auto locSavedIndex = locBunchIndexMap.find(locBeamBunches)->second;
	if(dDebugLevel >= 20)
	{
		cout << "Get_ResumeAtIndex_Combos: previous combo, bunches, saved index" << ", " << locPreviousCombo << ", ";
		for(auto& locBunch : locBeamBunches)
			cout << locBunch << ", ";
		cout << locSavedIndex << endl;
	}
	return locSavedIndex + 1;
}

inline bool DSourceComboer::Get_IsComboingZIndependent(const JObject* locObject, Particle_t locPID) const
{
	//make this virtual!
	if(ParticleCharge(locPID) != 0)
		return true;

	//actually, everything is z-dependent for massive neutrals.
	//however the momentum is SO z-dependent, that we can't cut on it until the end when we have the final vertex, AFTER comboing
	//a mere z-bin is not enough.
	//So, as far as COMBOING is concerned, massive neutrals are Z-INDEPENDENT
	if(ParticleMass(locPID) > 0.0)
		return true;

	auto locNeutralShower = static_cast<const DNeutralShower*>(locObject);
	return (locNeutralShower->dDetectorSystem == SYS_FCAL);
}

inline DSourceCombosByUse_Large& DSourceComboer::Get_CombosSoFar(ComboingStage_t locComboingStage, Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo)
{
	//THE INPUT locChargedCombo MUST BE:
	//If reading from: Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get
	//If saving to (you are making a mixed): Whatever charged combo you are about to combo horizontally with to make this new, mixed combo

	//NOTE: If on mixed stage, it is NOT valid to get fully-charged combos from here! In fact, what you want is probably the input combo!
	if(locComboingStage == d_ChargedStage)
		return dSourceCombosByUse_Charged;
	else if(locChargeContent_SearchForUse == d_Neutral)
		return dMixedCombosByUseByChargedCombo[nullptr]; //if fully neutral, then the charged combo doesn't matter: only matters for mixing charges
	return dMixedCombosByUseByChargedCombo[locChargedCombo];
}

inline DSourceCombosByBeamBunchByUse& DSourceComboer::Get_SourceCombosByBeamBunchByUse(Charge_t locChargeContent_SearchForUse, const DSourceCombo* locChargedCombo)
{
//shouldn't ever be called if charged
	//THE INPUT locChargedCombo MUST BE:
	//If reading from: Whatever charged combo you PREVIOUSLY comboed horizontally with to make the combos you're trying to get
	//If saving to (you are making a mixed): Whatever charged combo you just comboed horizontally with to make this new, mixed combo
	if(locChargeContent_SearchForUse == d_Neutral)
		return dSourceCombosByBeamBunchByUse[nullptr];
	return dSourceCombosByBeamBunchByUse[locChargedCombo];
}

inline bool DSourceComboer::Cut_dEdx(Particle_t locPID, DetectorSystem_t locSystem, double locP, double locdEdx)
{
	ddEdxValueMap[locPID][locSystem].emplace_back(locP, locdEdx);
	if(ddEdxCutMap[locPID].find(locSystem) == ddEdxCutMap[locPID].end())
		return true;

	auto locCutPair = ddEdxCutMap[locPID][locSystem];
//cout << "PID, p, dedx, cut value low/high, cut result: " << locPID << ", " << locP << ", " << locdEdx << ", " << locCutPair.first->Eval(locP) << ", " << locCutPair.second->Eval(locP) << ", " << ((locdEdx >= locCutPair.first->Eval(locP)) && (locdEdx <= locCutPair.second->Eval(locP))) << endl;
	return ((locdEdx >= locCutPair.first->Eval(locP)) && (locdEdx <= locCutPair.second->Eval(locP)));
}

inline bool DSourceComboer::Cut_EOverP(Particle_t locPID, DetectorSystem_t locSystem, double locP, double locEOverP)
{
	dEOverPValueMap[locPID][locSystem].emplace_back(locP, locEOverP);
	if(dEOverPCutMap[locPID].find(locSystem) == dEOverPCutMap[locPID].end())
		return true;

	auto locCutFunc = dEOverPCutMap[locPID][locSystem];
	return (IsLepton(locPID) == (locEOverP >= locCutFunc->Eval(locP)));
}

inline DSourceComboer::~DSourceComboer(void)
{
	//no need for a resource pool for these objects, as they will exist for the length of the program
	Fill_SurvivalHistograms();
	if(dDebugLevel >= 5)
	{
		Print_NumCombosByUse(); //for the final event

		cout << "FINAL Num combos by use (charged):" << endl;
		for(const auto& locNumCombosByUsePair : dNumMixedCombosMap_Charged)
		{
			cout << locNumCombosByUsePair.second << " of ";
			Print_SourceComboUse(locNumCombosByUsePair.first);
		}
		cout << "FINAL Num combos by use (neutral/mixed):" << endl;
		for(const auto& locNumCombosByUsePair : dNumMixedCombosMap_Mixed)
		{
			cout << locNumCombosByUsePair.second << " of ";
			Print_SourceComboUse(locNumCombosByUsePair.first);
		}

	}
	for(auto locComboInfo : dSourceComboInfos)
		delete locComboInfo;
	delete dSourceComboVertexer;
	delete dSourceComboP4Handler;
	delete dSourceComboTimeHandler;
	delete dParticleComboCreator;
}


/*********************************************************** INLINE NAMESPACE-SCOPE FUNCTION DEFINITIONS ************************************************************/

} //end DAnalysis namespace

#endif // DSourceComboer_h
