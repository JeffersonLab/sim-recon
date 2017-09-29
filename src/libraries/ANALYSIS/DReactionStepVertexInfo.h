#ifndef DReactionStepVertexInfo_h
#define DReactionStepVertexInfo_h

#include <set>
#include <map>
#include <iostream>
#include <vector>
#include <algorithm>

#include "particleType.h"
#include "ANALYSIS/DReaction.h"

using namespace std;

//reaction-based vertex info
namespace DAnalysis
{

enum DReactionState_t
{
	d_InitialState = 0,
	d_FinalState,
	d_EitherState
};

class DReactionStepVertexInfo
{
	public:

		//CONSTRUCTORS
		DReactionStepVertexInfo(void);
		void Set_Members(const DReaction* locReaction, size_t locStartStepIndex);
		void Reset(void);

		//ADDING STEPS & PARTICLES
		void Add_ReactionStep(size_t locStepIndex);
		void Set_ParticleIndices(bool locFitFlag, const vector<pair<int, int>>& locFullConstrainParticles, const vector<pair<int, int>>& locDecayingParticles,
				const vector<pair<int, int>>& locOnlyConstrainTimeParticles, const vector<pair<int, int>>& locNoConstrainParticles);

		//REGISTER DECAY VERTICES
		void Register_DecayingNoConstrainUseVertex(bool locFitFlag, const pair<int, int>& locDecayingNoConstrainPair, const DReactionStepVertexInfo* locVertexInfo);
		void Register_DecayingParticleConstraints(bool locFitFlag, const vector<pair<int, int>>& locNoConstrainDecayingParticles,
				const map<pair<int, int>, const DReactionStepVertexInfo*>& locFullConstrainDecayingParticles);
		void Set_DanglingVertexFlag(bool locIsDanglingVertexFlag){dIsDanglingVertexFlag = locIsDanglingVertexFlag;}
		void Set_FittableVertexFlag(bool locIsFittableVertexFlag){dIsFittableVertexFlag = locIsFittableVertexFlag;}
		void Set_IsInclusiveVertexFlag(bool locIsInclusiveVertexFlag){dIsInclusiveVertexFlag = locIsInclusiveVertexFlag;}

		//GET REACTION INFO
		const DReaction* Get_Reaction(void) const{return dReaction;}
		vector<size_t> Get_StepIndices(void) const{return dReactionStepIndices;}
		bool Get_ProductionVertexFlag(void) const{return dIsProductionVertexFlag;}
		bool Get_DanglingVertexFlag(void) const{return dIsDanglingVertexFlag;}
		bool Get_FittableVertexFlag(void) const{return dIsFittableVertexFlag;}
		bool Get_IsInclusiveVertexFlag(void) const{return dIsInclusiveVertexFlag;}

		//GET PARTICLES
		vector<pair<int, int>> Get_FullConstrainParticles(bool locFitFlag, DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges, bool locIncludeDecayingFlag = true) const;
		vector<pair<int, int>> Get_DecayingParticles(DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges) const;
		vector<pair<int, int>> Get_OnlyConstrainTimeParticles(void) const{return dOnlyConstrainTimeParticles;}
		vector<pair<int, int>> Get_NoConstrainParticles(bool locFitFlag, DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges,
				bool locIncludeDecayingFlag = true, bool locIncludeMissingFlag = true, bool locIncludeTargetFlag = true) const;
		vector<pair<int, int>> Get_MissingParticles(DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges) const;

		//sorted decaying
		map<pair<int, int>, const DReactionStepVertexInfo*> Get_DecayingParticles_NoConstrain(bool locFitFlag) const{return dDecayingParticles_NoConstrain.find(locFitFlag)->second;}
		map<pair<int, int>, const DReactionStepVertexInfo*> Get_DecayingParticles_FullConstrain(bool locFitFlag) const{return dDecayingParticles_FullConstrain.find(locFitFlag)->second;}

		//independent of state
		vector<pair<int, int>> Get_Particles(DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges,
				bool locIncludeDecayingFlag = true, bool locIncludeMissingFlag = true, bool locIncludeTargetFlag = true) const;

		//parent vertex info
		void Set_ParentVertexInfo(const DReactionStepVertexInfo* locStepVertexInfo){dParentVertexInfo = locStepVertexInfo;}
		const DReactionStepVertexInfo* Get_ParentVertexInfo(void) const{return dParentVertexInfo;}

	private:

		//MAP BOOL: KinFitFlag: true for kinfit particle information, false for reconstruction particle information
		//The difference:
		//For the production vertex, if the beamline is not used in the fit (due to errors = 0), and there's only one charged track:
		//The vertex is reconstructable, but it's not fittable

		//PRIVATE METHODS
		vector<pair<int, int>> Filter_Particles(vector<pair<int, int>> locParticles, DReactionState_t locState, Charge_t locCharge,
				bool locIncludeDecayingFlag = true, bool locIncludeMissingFlag = true, bool locIncludeTargetFlag = true) const;

		//REACTION SUMMARY INFO
		const DReaction* dReaction = nullptr; //only the first reaction is stored, though there may be several represented by this object (identical channel content)
		vector<size_t> dReactionStepIndices; //in order from smallest to largest
		bool dIsProductionVertexFlag = false;
		bool dIsInclusiveVertexFlag = false;

		//PARTICLE INFO //sorted!
		//pair: step, particle indices (including all: beam, target, decaying, detected, and missing)
		map<bool, vector<pair<int, int>>> dFullConstrainParticles; //detected charged tracks & beam, decaying when registered
		vector<pair<int, int>> dOnlyConstrainTimeParticles; //detected photons
		map<bool, vector<pair<int, int>>> dNoConstrainParticles; //missing, massive neutrals, decaying when registered, target

		//DECAY INFO
		//Note, decaying particles that decay in-place at this vertex (e.g. pi0) will only appear once: with their "final-state" indices
		//If the decaying particle has a detached vertex, then the indices reported here are the ones where it appears for the vertex (initial/final state)
		vector<pair<int, int>> dDecayingParticles; //all, whether used to constrain or not
		map<bool, map<pair<int, int>, const DReactionStepVertexInfo*>> dDecayingParticles_NoConstrain; //vertex-info: where it is used to constrain (nullptr if not)
		map<bool, map<pair<int, int>, const DReactionStepVertexInfo*>> dDecayingParticles_FullConstrain; //vertex-info: where it was defined

		//DANGLING
		//is it dangling? dangling = vertex indeterminable, even if not fitting & with all particle information
		//if is true, then vertex parent is either:
			//in dDecayingParticles_NoConstrain if it's not empty (at most one will have non-null info), or is center of target
		bool dIsDanglingVertexFlag = false;
		bool dIsFittableVertexFlag = true;
		const DReactionStepVertexInfo* dParentVertexInfo = nullptr; //null if production vertex
};

/****************************************************** NAMESPACE-SCOPE NON-INLINE FUNCTION DECLARATIONS *******************************************************/

string Build_VertexConstraintString(const DReactionStepVertexInfo* locVertexInfo, bool locSpacetimeFitFlag);
void Print_ReactionStepVertexInfo(const DReactionStepVertexInfo* locStepInfo);

/************************************************************** CONSTRUCTORS & OPERATORS ***************************************************************/

inline DReactionStepVertexInfo::DReactionStepVertexInfo(void)
{
	dFullConstrainParticles.emplace(true, vector<pair<int, int>>{});
	dFullConstrainParticles.emplace(false, vector<pair<int, int>>{});
	dNoConstrainParticles.emplace(true, vector<pair<int, int>>{});
	dNoConstrainParticles.emplace(false, vector<pair<int, int>>{});

	dDecayingParticles_NoConstrain.emplace(true, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_NoConstrain.emplace(false, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_FullConstrain.emplace(true, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_FullConstrain.emplace(false, map<pair<int, int>, const DReactionStepVertexInfo*>{});
}

inline void DReactionStepVertexInfo::Set_Members(const DReaction* locReaction, size_t locStartStepIndex)
{
	dReaction = locReaction;
	dReactionStepIndices = {locStartStepIndex};
	dIsProductionVertexFlag = ((locStartStepIndex == 0) && DAnalysis::Get_IsFirstStepBeam(locReaction));
}

inline void DReactionStepVertexInfo::Reset(void)
{
	//REACTION SUMMARY INFO
	dReaction = nullptr;
	dReactionStepIndices.clear();
	dIsProductionVertexFlag = false;
	dIsInclusiveVertexFlag = false;

	//PARTICLE INFO
	dFullConstrainParticles.clear();
	dOnlyConstrainTimeParticles.clear();
	dNoConstrainParticles.clear();

	//DECAY INFO
	dDecayingParticles.clear();
	dDecayingParticles_NoConstrain.clear();
	dDecayingParticles_FullConstrain.clear();

	//EMPLACE
	dFullConstrainParticles.emplace(true, vector<pair<int, int>>{});
	dFullConstrainParticles.emplace(false, vector<pair<int, int>>{});
	dNoConstrainParticles.emplace(true, vector<pair<int, int>>{});
	dNoConstrainParticles.emplace(false, vector<pair<int, int>>{});

	dDecayingParticles_NoConstrain.emplace(true, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_NoConstrain.emplace(false, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_FullConstrain.emplace(true, map<pair<int, int>, const DReactionStepVertexInfo*>{});
	dDecayingParticles_FullConstrain.emplace(false, map<pair<int, int>, const DReactionStepVertexInfo*>{});

	//FLAGS
	dIsDanglingVertexFlag = false;
	dIsFittableVertexFlag = true;
	dParentVertexInfo = nullptr;
}

/************************************************************** INLINE FUNCTIONS ***************************************************************/

inline void DReactionStepVertexInfo::Add_ReactionStep(size_t locStepIndex)
{
	dReactionStepIndices.push_back(locStepIndex);
	std::sort(dReactionStepIndices.begin(), dReactionStepIndices.end()); //just in case
}

inline void DReactionStepVertexInfo::Register_DecayingNoConstrainUseVertex(bool locFitFlag, const pair<int, int>& locDecayingNoConstrainPair, const DReactionStepVertexInfo* locVertexInfo)
{
	dDecayingParticles_NoConstrain[locFitFlag][locDecayingNoConstrainPair] = locVertexInfo;
}

inline vector<pair<int, int>> DReactionStepVertexInfo::Get_FullConstrainParticles(bool locFitFlag, DReactionState_t locState, Charge_t locCharge, bool locIncludeDecayingFlag) const
{
	return Filter_Particles(dFullConstrainParticles.find(locFitFlag)->second, locState, locCharge, locIncludeDecayingFlag);
}

inline vector<pair<int, int>> DReactionStepVertexInfo::Get_DecayingParticles(DReactionState_t locState, Charge_t locCharge) const
{
	return Filter_Particles(dDecayingParticles, locState, locCharge);
}

inline vector<pair<int, int>> DReactionStepVertexInfo::Get_NoConstrainParticles(bool locFitFlag, DReactionState_t locState, Charge_t locCharge,
		bool locIncludeDecayingFlag, bool locIncludeMissingFlag, bool locIncludeTargetFlag) const
{
	return Filter_Particles(dNoConstrainParticles.find(locFitFlag)->second, locState, locCharge, locIncludeDecayingFlag, locIncludeMissingFlag, locIncludeTargetFlag);
}

inline vector<pair<int, int>> DReactionStepVertexInfo::Get_MissingParticles(DReactionState_t locState, Charge_t locCharge) const
{
	vector<pair<int, int>> locParticles = Filter_Particles(dNoConstrainParticles.find(false)->second, locState, locCharge, false, true, false);

	auto Check_NotMissing = [&](const pair<int, int>& locIndices) -> bool
	{return (locIndices.second != dReaction->Get_ReactionStep(locIndices.first)->Get_MissingParticleIndex());};

	locParticles.erase(std::remove_if(locParticles.begin(), locParticles.end(), Check_NotMissing), locParticles.end());
	return locParticles;
}

inline vector<pair<int, int>> DReactionStepVertexInfo::Get_Particles(DReactionState_t locState, Charge_t locCharge, bool locIncludeDecayingFlag, bool locIncludeMissingFlag, bool locIncludeTargetFlag) const
{
	//ASSUMES Object has been fully created before calling
	vector<pair<int, int>> locParticles = dFullConstrainParticles.find(true)->second;
	locParticles.insert(locParticles.end(), dOnlyConstrainTimeParticles.begin(), dOnlyConstrainTimeParticles.end());
	locParticles.insert(locParticles.end(), dNoConstrainParticles.find(true)->second.begin(), dNoConstrainParticles.find(true)->second.end());
	std::sort(locParticles.begin(), locParticles.end());
	return Filter_Particles(locParticles, locState, locCharge, locIncludeDecayingFlag, locIncludeMissingFlag, locIncludeTargetFlag);
}

} //end DAnalysis namespace

#endif // DReactionStepVertexInfo_h
