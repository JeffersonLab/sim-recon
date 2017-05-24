#ifndef DParticleComboCreator_h
#define DParticleComboCreator_h

using namespace std;
using namespace jana;

#include "JANA/JEventLoop.h"
#include "PID/DEventRFBunch.h"
#include "PID/DChargedTrackHypothesis_factory.h"
#include "PID/DNeutralParticleHypothesis_factory.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboStep.h"

#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DSourceComboVertexer.h"

namespace DAnalysis
{

class DParticleComboCreator
{
	public:
		void DParticleComboCreator(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboTimeHandler* locSourceComboTimeHandler, const DSourceComboVertexer* dSourceComboVertexer);

		const DParticleCombo* Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift);

	private:

		//UTILITIES
		const DSourceComboer* dSourceComboer;
		const DSourceComboTimeHandler* dSourceComboTimeHandler;
		const DSourceComboVertexer* dSourceComboVertexer;

		//FACTORIES
		DNeutralParticleHypothesis_factory* dNeutralParticleHypothesisFactory;
		DChargedTrackHypothesis_factory* dChargedTrackHypothesisFactory;

		//CREATED OBJECT MAPS
		unordered_map<tuple<const DReactionStep*, const DSourceCombo*, const DKinematicData*>, const DParticleComboStep*> dComboStepMap; //kindata is beam (null if not in step): for vertex
		map<int, const DEventRFBunch*> dRFBunchMap;

		//RESOURCE POOLS
		DResourcePool<DEventRFBunch> dResourcePool_EventRFBunch;
		DResourcePool<DParticleCombo> dResourcePool_ParticleCombo;
		DResourcePool<DParticleComboStep> dResourcePool_ParticleComboStep;
};

}

#endif // DParticleComboCreator_h
