#ifndef DParticleComboCreator_h
#define DParticleComboCreator_h

using namespace std;
using namespace jana;

#include <unordered_map>

#include "JANA/JEventLoop.h"
#include "PID/DEventRFBunch.h"
#include "PID/DChargedTrackHypothesis_factory.h"
#include "PID/DNeutralParticleHypothesis_factory.h"
#include "PID/DBeamPhoton_factory.h"
#include "PID/DParticleID.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboStep.h"

#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DSourceComboVertexer.h"

namespace DAnalysis
{

class DSourceComboer;
class DParticleComboCreator
{
	public:
		void DParticleComboCreator(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboTimeHandler* locSourceComboTimeHandler, const DSourceComboVertexer* dSourceComboVertexer);

		const DParticleCombo* Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift, DKinFitType locKinFitType);

		void Reset(void);

	private:

		bool Get_CreateNeutralErrorMatrixFlag_Combo(const DReactionVertexInfo* locReactionVertexInfo, DKinFitType locKinFitType);

		//CREATE PARTICLES
		const DChargedTrackHypothesis* Create_ChargedHypo(const DChargedTrack* locChargedTrack, Particle_t locPID, double locPropagatedRFTime, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryFullCombo, const DKinematicData* locBeamParticle);
		const DBeamPhoton* Build_BeamPhoton_KinFit(const DBeamPhoton* locBeamPhoton, DKinFitParticle* locKinFitParticle);

		TMatrixFSym dVertexCovMatrix;
		unordered_map<const DReactionVertexInfo*, bool> dDanglingNeutralsFlagMap;

		//UTILITIES
		const DSourceComboer* dSourceComboer;
		const DSourceComboTimeHandler* dSourceComboTimeHandler;
		const DSourceComboVertexer* dSourceComboVertexer;
		const DParticleID* dParticleID;

		//FACTORIES
		DNeutralParticleHypothesis_factory* dNeutralParticleHypothesisFactory;
		DChargedTrackHypothesis_factory* dChargedTrackHypothesisFactory;
		DBeamPhoton_factory* dBeamPhotonfactory;

		//CREATED OBJECT MAPS
		unordered_map<tuple<const DSourceCombo*, bool, bool, const DSourceCombo*, const DKinematicData*>, const DParticleComboStep*> dComboStepMap; //kindata is beam (null if not in step): for vertex
		unordered_map<int, const DEventRFBunch*> dRFBunchMap;
		unordered_map<tuple<const DChargedTrack*, Particle_t, int, bool, bool, const DSourceCombo*, const DSourceCombo*, const DKinematicData*>, DChargedTrackHypothesis*> dChargedHypoMap;
		unordered_map<tuple<const DNeutralShower*, Particle_t, int, bool, const DSourceCombo*, const DSourceCombo*, const DKinematicData*>, DNeutralParticleHypothesis*> dNeutralHypoMap;
		unordered_map<tuple<const DReactionVertexInfo*, const DSourceCombo*, const DKinematicData*, int, bool>, DParticleCombo*> dComboMap;

		//RESOURCE POOLS
		DResourcePool<DEventRFBunch> dResourcePool_EventRFBunch;
		DResourcePool<DParticleCombo> dResourcePool_ParticleCombo;
		DResourcePool<DParticleComboStep> dResourcePool_ParticleComboStep;
};

}

#endif // DParticleComboCreator_h
