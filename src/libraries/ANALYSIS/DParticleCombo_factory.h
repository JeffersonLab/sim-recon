#ifndef _DParticleCombo_factory_
#define _DParticleCombo_factory_

#include <deque>
#include <vector>

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "particleType.h"
#include "TMatrixDSym.h"

#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DKinematicData.h"
#include "PID/DBeamPhoton.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboBlueprint.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "KINFITTER/DKinFitParticle.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "ANALYSIS/DAnalysisResults.h"

using namespace std;
using namespace jana;

class DParticleCombo_factory : public jana::JFactory<DParticleCombo>
{
	public:
		DParticleCombo_factory(){};
		~DParticleCombo_factory(){};

		void Reset_Pools(void);

		size_t Get_ParticleComboStepPoolSize(void) const{return dParticleComboStepPool_All.size();};
		size_t Get_KinematicDataPoolSize(void) const{return dKinematicDataPool_All.size();};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop* locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DParticleCombo* Check_IsDuplicateCombo(const map<const DParticleCombo*, DParticleCombo*>& locNewParticleComboMap, const DParticleCombo* locParticleCombo);

		void Set_DecayingParticles(const DParticleCombo* locNewParticleCombo, const DParticleCombo* locOldParticleCombo, size_t locStepIndex, DParticleComboStep* locNewParticleComboStep, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		DKinFitParticle* Get_DecayingParticle(const DParticleCombo* locOldParticleCombo, size_t locComboStepIndex, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		bool Search_ForParticleInDecay(const DKinFitChain* locKinFitChain, size_t locStepToSearch, DKinFitParticle* locParticleToFind);

		const DChargedTrackHypothesis* Get_ChargedHypothesis(const DParticleCombo* locParticleCombo, const vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses, const DKinematicData* locKinematicData_Measured) const;
		const DNeutralParticleHypothesis* Get_NeutralHypothesis(const DParticleCombo* locParticleCombo, const vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses, const DKinematicData* locKinematicData_Measured) const;

		void Set_SpacetimeVertex(const DParticleCombo* locNewParticleCombo, DParticleComboStep* locNewParticleComboStep, size_t locStepIndex, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain) const;
		void Reset_Data(void);

		DKinematicData* Build_KinematicData(Particle_t locPID, DKinFitParticle* locKinFitParticle);

		DKinematicData* Get_KinematicDataResource(void);
		DParticleComboStep* Get_ParticleComboStepResource(void);

		DKinFitUtils_GlueX* dKinFitUtils;

		vector<DParticleCombo*> dCreatedParticleCombos;

		deque<DParticleComboStep*> dParticleComboStepPool_All;
		deque<DParticleComboStep*> dParticleComboStepPool_Available;

		deque<DKinematicData*> dKinematicDataPool_All;
		deque<DKinematicData*> dKinematicDataPool_Available;

		size_t MAX_DParticleComboStepPoolSize;
		size_t MAX_DKinematicDataPoolSize;
};

#endif // _DParticleCombo_factory_


