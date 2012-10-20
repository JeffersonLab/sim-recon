// $Id$
//
//    File: DParticleCombo_factory_PreKinFit.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleCombo_factory_PreKinFit_
#define _DParticleCombo_factory_PreKinFit_

#include <deque>
#include <vector>
#include <map>

#include "JANA/JFactory.h"
#include "JANA/JObject.h"

#include "particleType.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DKinematicData.h"
#include "PID/DChargedTrack.h"
#include "PID/DBeamPhoton.h"
#include "PID/DNeutralShower.h"

#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DParticleComboBlueprint.h"

class DParticleCombo_factory_PreKinFit : public jana::JFactory<DParticleCombo>
{
	public:
		DParticleCombo_factory_PreKinFit(){};
		~DParticleCombo_factory_PreKinFit(){};
		const char* Tag(void){return "PreKinFit";}

		void Reset_Pools(void);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DKinematicData* Get_DetectedParticle(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locParticleIndex, vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses_Reaction, vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses, const JObject*& locSourceObject);
		DKinematicData* Create_Target(Particle_t locPID);
		DBeamPhoton* Create_BeamPhoton(void); //for MC only!

		double dMaxPhotonRFTimeDifference;

		DParticleComboStep* Clone_ParticleComboStep(const DParticleComboStep* locParticleComboStep);

		void Reset_KinematicData(DKinematicData* locKinematicData);

		DParticleComboStep* Get_ParticleComboStepResource(void);
		DKinematicData* Get_KinematicDataResource(void);
		DBeamPhoton* Get_BeamPhotonResource(void);

		deque<DParticleComboStep*> dParticleComboStepPool_All;
		deque<DParticleComboStep*> dParticleComboStepPool_Available;

		deque<DKinematicData*> dKinematicDataPool_All;
		deque<DKinematicData*> dKinematicDataPool_Available;

		deque<DBeamPhoton*> dBeamPhotonPool_All; //for MC only!
		deque<DBeamPhoton*> dBeamPhotonPool_Available;

		map<const DParticleComboBlueprintStep*, const DParticleComboStep*> dComboBlueprintStepMap;

		size_t MAX_DParticleComboStepPoolSize;
		size_t MAX_DKinematicDataPoolSize;
		size_t MAX_DBeamPhotonPoolSize;

		bool dVertexZCutFlag;
		double dMinVertexZ;
		double dMaxVertexZ;

		double dMinChargedPIDFOM;
		double dMaxTrackingChiSqPerDF;
};

#endif // _DParticleCombo_factory_PreKinFit_

