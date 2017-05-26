#ifndef DReactionVertexInfo_h
#define DReactionVertexInfo_h

#include <memory>
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
		DReactionVertexInfo(const DReaction* locReaction, const vector<shared_ptr<DReactionStepVertexInfo>>& locStepVertexInfos);

		//SETTERS
		void Add_Reaction(const DReaction* locReaction){return dReactions.insert(locReaction);}

		//GETTERS
		const DReaction* Get_Reaction(void) const{return *dReactions.begin();} //since their channels are identical, any one will do (if used correctly)
		vector<const DReaction*> Get_Reactions(void) const{return dReactions;}
		vector<shared_ptr<const DReactionStepVertexInfo>> Get_StepVertexInfos(void) const{return dStepVertexInfos;}
		shared_ptr<const DReactionStepVertexInfo> Get_StepVertexInfo(size_t locStepIndex) const{return dVertexInfoMap.at(locStepIndex);}

	private:

		//these all have identical channel content: particles & steps must be in the same order, although actions, etc. may be different
		unordered_set<const DReaction*> dReactions;

		vector<shared_ptr<const DReactionStepVertexInfo>> dStepVertexInfos; //in order of construction dependency
		unordered_map<size_t, shared_ptr<const DReactionStepVertexInfo>> dVertexInfoMap; //key is step index
};

inline DReactionVertexInfo::DReactionVertexInfo(const DReaction* locReaction, const vector<shared_ptr<DReactionStepVertexInfo>>& locStepVertexInfos) :
		dReactions({locReaction})
{
	//transform into vector of shared_ptr containing const pointers
	auto locConstify = [](shared_ptr<DReactionStepVertexInfo>& locVertexInfo) {return std::const_pointer_cast<const DReactionStepVertexInfo>(locVertexInfo);};
	std::transform(locStepVertexInfos.begin(), locStepVertexInfos.end(), std::back_inserter(dStepVertexInfos), locConstify);

	//build the step index map
	for(auto locVertexInfo : dStepVertexInfos)
	{
		for(auto locStepIndex : locVertexInfo->Get_StepIndices())
			dVertexInfoMap.emplace(locStepIndex, locVertexInfo);
	}
}

//NAMESPACE SCOPE FUNCTIONS
inline vector<shared_ptr<const DReactionStepVertexInfo>> Get_StepVertexInfos_ReverseOrderByStep(const DReactionVertexInfo* locReactionVertexInfo)
{
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();

	//sort vertex infos in reverse-step order
	auto Comparator_ReverseOrderByStep = [](const DReactionStepVertexInfo* lhs, const DReactionStepVertexInfo* rhs) -> bool
		{return lhs->Get_StepIndices().front() > rhs->Get_StepIndices().front();}; // >: reverse order

	std::sort(locStepVertexInfos.begin(), locStepVertexInfos.end(), Comparator_ReverseOrderByStep);
	return locStepVertexInfos;
}

vector<pair<int, int>> Get_FullConstrainParticles(const DReactionVertexInfo* locReactionVertexInfo, DReactionState_t locState = d_EitherState, Charge_t locCharge = d_AllCharges, bool locIncludeDecayingFlag = true)
{
	//no longer sorted!
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	vector<pair<int, int>> locAllFullConstrainParticles;
	for(auto& locStepVertexInfo : locStepVertexInfos)
	{
		auto locFullConstrainParticles = locStepVertexInfo->Get_FullConstrainParticles(locState, locCharge, locIncludeDecayingFlag);
		locAllFullConstrainParticles.insert(locAllFullConstrainParticles.end(), locFullConstrainParticles.begin(), locFullConstrainParticles.end());
	}
	return locAllFullConstrainParticles;
}

vector<pair<int, int>> Get_OnlyConstrainTimeParticles(const DReactionVertexInfo* locReactionVertexInfo)
{
	//no longer sorted!
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	vector<pair<int, int>> locAllOnlyConstrainTimeParticles;
	for(auto& locStepVertexInfo : locStepVertexInfos)
	{
		auto locOnlyConstrainTimeParticles = locStepVertexInfo->Get_OnlyConstrainTimeParticles();
		locAllOnlyConstrainTimeParticles.insert(locAllFullConstrainParticles.end(), locOnlyConstrainTimeParticles.begin(), locOnlyConstrainTimeParticles.end());
	}
	return locAllOnlyConstrainTimeParticles;
}

} //end DAnalysis namespace

#endif // DReactionVertexInfo_h
