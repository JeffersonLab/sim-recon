#ifndef _DParticleCombo_factory_Thrown_
#define _DParticleCombo_factory_Thrown_

#include <iostream>
#include <deque>

#include "JANA/JFactory.h"
#include "particleType.h"

#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DParticleComboBlueprintStep.h"
#include "PID/DMCReaction.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"

class DAnalysisUtilities;

using namespace jana;
using namespace std;

class DParticleCombo_factory_Thrown : public jana::JFactory<DParticleCombo>
{
	public:
		DParticleCombo_factory_Thrown(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DParticleCombo_factory_Thrown(){};
		const char* Tag(void){return "Thrown";}

		DParticleCombo* Build_ThrownCombo(JEventLoop* locEventLoop, const DReaction* locThrownReaction, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DAnalysisUtilities* dAnalysisUtilities;

		DParticleComboStep* Get_ParticleComboStepResource(void);
		DParticleComboBlueprintStep* Get_ParticleComboBlueprintStepResource(void);

		deque<DParticleComboStep*> dParticleComboStepPool_All;
		deque<DParticleComboStep*> dParticleComboStepPool_Available;

		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_All;
		deque<DParticleComboBlueprintStep*> dParticleComboBlueprintStepPool_Available;

		size_t MAX_dParticleComboStepPoolSize;
};

#endif // _DParticleCombo_factory_Thrown_

