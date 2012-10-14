#ifndef _DParticleCombo_factory_
#define _DParticleCombo_factory_

#include <deque>
#include <vector>

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "TMatrixDSym.h"

#include "PID/DKinFit.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboBlueprint.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "ANALYSIS/DKinFitParticle.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DKinFitter_GlueX.h"
#include "ANALYSIS/DAnalysisResults.h"

using namespace std;
using namespace jana;

class DParticleCombo_factory : public jana::JFactory<DParticleCombo>
{
	public:
		DParticleCombo_factory(){};
		~DParticleCombo_factory(){};

		void Reset_Pools(void);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DKinFitter_GlueX dKinFitter;

		DParticleCombo* Clone_ParticleCombo(const DParticleCombo* locParticleCombo);
		DKinematicData* Build_KinematicData(Particle_t locPID, const DKinFitParticle* locKinFitParticle);

		DKinematicData* Get_KinematicDataResource(void);
		DParticleComboStep* Get_ParticleComboStepResource(void);

		void Reset_KinematicData(DKinematicData* locKinematicData);

		deque<DParticleComboStep*> dParticleComboStepPool_All;
		deque<DParticleComboStep*> dParticleComboStepPool_Available;

		deque<DKinematicData*> dKinematicDataPool_All;
		deque<DKinematicData*> dKinematicDataPool_Available;

		size_t MAX_DParticleComboStepPoolSize;
		size_t MAX_DKinematicDataPoolSize;
};

#endif // _DParticleCombo_factory_

