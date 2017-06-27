#ifndef _DNeutralParticle_factory_Combo_
#define _DNeutralParticle_factory_Combo_

#include <string>
#include <vector>
#include <set>

#include <JANA/JFactory.h>
#include "PID/DEventRFBunch.h"
#include <PID/DNeutralParticle.h>
#include "PID/DVertex.h"
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DNeutralShower.h>
#include "PID/DNeutralParticleHypothesis_factory.h"
#include "ANALYSIS/DReaction.h"

using namespace std;
using namespace jana;

class DNeutralParticle_factory_Combo : public jana::JFactory<DNeutralParticle>
{
	public:
		DNeutralParticle_factory_Combo(){};
		~DNeutralParticle_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.

		string dShowerSelectionTag;
		set<Particle_t> dNeutralPIDs;
		vector<DNeutralParticleHypothesis*> dCreatedHypotheses;
		DNeutralParticleHypothesis_factory* dNeutralParticleHypothesisFactory;
};

#endif // _DNeutralParticle_factory_Combo_
