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

	//GETTERS
	const DReaction* Get_Reaction(void) const{return dReaction;}
	vector<shared_ptr<const DReactionStepVertexInfo>> Get_StepVertexInfos(void) const{return dStepVertexInfos;}
	shared_ptr<const DReactionStepVertexInfo> Get_StepVertexInfo(size_t locStepIndex) const{return dVertexInfoMap.at(locStepIndex);}

private:
	const DReaction* dReaction;
	vector<shared_ptr<const DReactionStepVertexInfo>> dStepVertexInfos; //in order of construction dependency
	unordered_map<size_t, shared_ptr<const DReactionStepVertexInfo>> dVertexInfoMap; //key is step index
};

inline DReactionVertexInfo::DReactionVertexInfo(const DReaction* locReaction, const vector<shared_ptr<DReactionStepVertexInfo>>& locStepVertexInfos) :
		dReaction(locReaction)
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

} //end DAnalysis namespace

#endif // DReactionVertexInfo_h
