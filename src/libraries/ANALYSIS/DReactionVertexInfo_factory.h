#ifndef DReactionVertexInfo_factory_h
#define DReactionVertexInfo_factory_h

#include <unordered_map>
#include <map>
#include <vector>
#include <algorithm>

#include "JANA/JFactory.h"

#include "particleType.h"
#include "DResourcePool.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DReactionVertexInfo.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"

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
		jerror_t fini(void)
		{
			for(auto locInfo : _data)
				delete locInfo;
			_data.clear();
			delete dResourcePool_ReactionStepVertexInfo;
			return NOERROR;
		}

		DReactionVertexInfo* Build_VertexInfo(const DReaction* locReaction);

		//SETUP
		DReactionStepVertexInfo* Setup_VertexInfo(const DReaction* locReaction, size_t locStepIndex, DReactionStepVertexInfo* locVertexInfo);

		//GROUPING
		void Group_VertexParticles(DReactionStepVertexInfo* locVertexInfo);
		vector<DReactionStepVertexInfo*> Link_Vertices(const DReaction* locReaction, vector<DReactionStepVertexInfo*> locVertexInfos, bool locFitFlag) const;
		bool Associate_DecayingParticles(bool locFitFlag, bool locLinkingFlag, DReactionStepVertexInfo* locVertexInfo, map<pair<int, int>, DReactionStepVertexInfo*>& locDefinedDecayingParticles) const;

		//not all reactions are stored here, just the first ones
		unordered_map<const DReaction*, DReactionVertexInfo*> dVertexInfoMap;

		//RESOURCE POOL
		DResourcePool<DReactionStepVertexInfo>* dResourcePool_ReactionStepVertexInfo = nullptr;

		DKinFitUtils_GlueX* dKinFitUtils;
};

} //end DAnalysis namespace

#endif // DReactionVertexInfo_factory_h
