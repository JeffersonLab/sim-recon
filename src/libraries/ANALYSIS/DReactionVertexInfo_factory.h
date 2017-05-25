#ifndef DReactionVertexInfo_factory_h
#define DReactionVertexInfo_factory_h

#include <memory>
#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>

#include "JANA/JFactory.h"

#include <particleType.h>
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionVertexInfo.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

class DReactionVertexInfo_factory : public jana::JFactory<DReactionVertexInfo>
{
	private:

		//PRIMARY FUNCTIONS
		jerror_t init(void);
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t locEventNumber);
		DReactionVertexInfo* Build_VertexInfo(const DReaction* locReaction) const;

		//SETUP
		shared_ptr<DReactionStepVertexInfo> Setup_VertexInfo(const DReaction* locReaction, size_t locStepIndex, DReactionStepVertexInfo* locVertexInfo) const;

		//GROUPING
		void Group_VertexParticles(DReactionStepVertexInfo* locVertexInfo);
		vector<shared_ptr<DReactionStepVertexInfo>> Link_Vertices(const DReaction* locReaction, vector<shared_ptr<DReactionStepVertexInfo>> locVertexInfos) const;
		bool Associate_DecayingParticles(bool locLinkingFlag, shared_ptr<DReactionStepVertexInfo>& locVertexInfo, map<pair<int, int>, shared_ptr<DReactionStepVertexInfo>>& locDefinedDecayingParticles) const;

		//not all reactions are stored here, just the first ones
		unordered_map<const DReaction*, DReactionVertexInfo*> dVertexInfoMap;
};

} //end DAnalysis namespace

#endif // DReactionVertexInfo_factory_h
