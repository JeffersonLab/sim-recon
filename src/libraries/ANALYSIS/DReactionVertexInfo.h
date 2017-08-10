#ifndef DReactionVertexInfo_h
#define DReactionVertexInfo_h

#include <unordered_map>
#include <vector>
#include <algorithm>

#include "JANA/JObject.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionStepVertexInfo.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

class DReactionVertexInfo : public JObject
{
	public:
		JOBJECT_PUBLIC(DReactionVertexInfo);

		//CONSTRUCTORS
		DReactionVertexInfo(void) = delete;
		DReactionVertexInfo(const DReaction* locReaction, const vector<DReactionStepVertexInfo*>& locStepVertexInfos);

		//SETTERS
		void Add_Reaction(const DReaction* locReaction){dReactions.push_back(locReaction);}

		//GETTERS
		const DReaction* Get_Reaction(void) const{return *dReactions.begin();} //since their channels are identical, any one will do (if used correctly)
		vector<const DReaction*> Get_Reactions(void) const{return dReactions;}
		vector<const DReactionStepVertexInfo*> Get_StepVertexInfos(void) const{return dStepVertexInfos;}
		vector<const DReactionStepVertexInfo*> Get_StepVertexInfos_StepOrder(void) const;
		const DReactionStepVertexInfo* Get_StepVertexInfo(size_t locStepIndex) const{return dVertexInfoMap.at(locStepIndex);}

	private:

		//these all have identical channel content: particles & steps must be in the same order, although actions, etc. may be different
		vector<const DReaction*> dReactions;

		vector<const DReactionStepVertexInfo*> dStepVertexInfos; //in order of construction dependency
		unordered_map<size_t, const DReactionStepVertexInfo*> dVertexInfoMap; //key is step index
};

inline void Print_ReactionVertexInfo(const DReactionVertexInfo* locReactionInfo)
{
	cout << "Reaction name: " << locReactionInfo->Get_Reaction()->Get_ReactionName() << endl;
	for(auto& locStepInfo : locReactionInfo->Get_StepVertexInfos())
		DAnalysis::Print_ReactionStepVertexInfo(locStepInfo);
}

inline DReactionVertexInfo::DReactionVertexInfo(const DReaction* locReaction, const vector<DReactionStepVertexInfo*>& locStepVertexInfos) :
		dReactions({locReaction})
{
	//transform into vector containing const pointers
	auto locConstify = [](DReactionStepVertexInfo* locVertexInfo) -> const DReactionStepVertexInfo* {return const_cast<const DReactionStepVertexInfo*>(locVertexInfo);};
	std::transform(locStepVertexInfos.begin(), locStepVertexInfos.end(), std::back_inserter(dStepVertexInfos), locConstify);

	//build the step index map
	for(auto locVertexInfo : dStepVertexInfos)
	{
		for(auto locStepIndex : locVertexInfo->Get_StepIndices())
			dVertexInfoMap.emplace(locStepIndex, locVertexInfo);
	}
}

inline vector<const DReactionStepVertexInfo*> DReactionVertexInfo::Get_StepVertexInfos_StepOrder(void) const
{
	vector<const DReactionStepVertexInfo*> locVertexInfos;
	for(const auto& locStepPair : dVertexInfoMap)
	{
		if(locStepPair.first == 0)
			locVertexInfos.push_back(locStepPair.second);
		else if(locStepPair.second != locVertexInfos.back())
			locVertexInfos.push_back(locStepPair.second);
	}
	return locVertexInfos;
}

//NAMESPACE SCOPE FUNCTIONS
inline vector<const DReactionStepVertexInfo*> Get_StepVertexInfos_ReverseOrderByStep(const DReactionVertexInfo* locReactionVertexInfo)
{
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();

	//sort vertex infos in reverse-step order
	auto Comparator_ReverseOrderByStep = [](const DReactionStepVertexInfo* lhs, const DReactionStepVertexInfo* rhs) -> bool
		{return lhs->Get_StepIndices().front() > rhs->Get_StepIndices().front();}; // >: reverse order

	std::sort(locStepVertexInfos.begin(), locStepVertexInfos.end(), Comparator_ReverseOrderByStep);
	return locStepVertexInfos;
}

inline vector<pair<int, int>> Get_FullConstrainParticles(const DReactionVertexInfo* locReactionVertexInfo, bool locFitFlag, DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges, bool locIncludeDecayingFlag = true)
{
	//no longer sorted!
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	vector<pair<int, int>> locAllFullConstrainParticles;
	for(auto& locStepVertexInfo : locStepVertexInfos)
	{
		auto locFullConstrainParticles = locStepVertexInfo->Get_FullConstrainParticles(locFitFlag, locState, locCharge, locIncludeDecayingFlag);
		locAllFullConstrainParticles.insert(locAllFullConstrainParticles.end(), locFullConstrainParticles.begin(), locFullConstrainParticles.end());
	}
	return locAllFullConstrainParticles;
}

inline vector<pair<int, int>> Get_OnlyConstrainTimeParticles(const DReactionVertexInfo* locReactionVertexInfo)
{
	//no longer sorted!
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	vector<pair<int, int>> locAllOnlyConstrainTimeParticles;
	for(auto& locStepVertexInfo : locStepVertexInfos)
	{
		auto locOnlyConstrainTimeParticles = locStepVertexInfo->Get_OnlyConstrainTimeParticles();
		locAllOnlyConstrainTimeParticles.insert(locAllOnlyConstrainTimeParticles.end(), locOnlyConstrainTimeParticles.begin(), locOnlyConstrainTimeParticles.end());
	}
	return locAllOnlyConstrainTimeParticles;
}

} //end DAnalysis namespace

#endif // DReactionVertexInfo_h
