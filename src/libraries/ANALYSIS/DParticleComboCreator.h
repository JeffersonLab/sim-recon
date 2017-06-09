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
		const DParticleCombo* Create_KinFitCombo_NewCombo(const DParticleCombo* locOrigCombo, const DReaction* locReaction, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain);

		void Reset(void);

	private:

		bool Get_CreateNeutralErrorMatrixFlag_Combo(const DReactionVertexInfo* locReactionVertexInfo, DKinFitType locKinFitType);

		//DECAYING PARTICLES, POST-KINFIT
		void Set_DecayingParticles(const DParticleCombo* locNewParticleCombo, const DParticleCombo* locOldParticleCombo, size_t locStepIndex, DParticleComboStep* locNewParticleComboStep, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		DKinFitParticle* Get_DecayingParticle(const DParticleCombo* locOldParticleCombo, size_t locComboStepIndex, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		bool Search_ForParticleInDecay(const DKinFitChain* locKinFitChain, size_t locStepToSearch, DKinFitParticle* locParticleToFind);

		//SPACETIME VERTEX, POST-KINFIT
		void Set_SpacetimeVertex(const DParticleCombo* locNewParticleCombo, DParticleComboStep* locNewParticleComboStep, size_t locStepIndex, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain) const;

		//CREATE PARTICLES
		const DChargedTrackHypothesis* Create_ChargedHypo(const DChargedTrack* locChargedTrack, Particle_t locPID, double locPropagatedRFTime, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryFullCombo, const DKinematicData* locBeamParticle);
		const DBeamPhoton* Create_BeamPhoton_KinFit(const DBeamPhoton* locBeamPhoton, const DKinFitParticle* locKinFitParticle);
		const DChargedTrackHypothesis* Create_ChargedHypo_KinFit(const DChargedTrackHypothesis* locOrigHypo, const DKinFitParticle* locKinFitParticle, DKinFitType locKinFitType);
		const DNeutralParticleHypothesis* Create_NeutralHypo_KinFit(const DNeutralParticleHypothesis* locOrigHypo, DKinFitParticle* locKinFitParticle, DKinFitType locKinFitType);
		DKinematicData* Build_KinematicData(DKinFitParticle* locKinFitParticle, DKinFitType locKinFitType, DVector3 locPreKinFitVertex);

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
		unordered_map<const DKinFitParticle*, DChargedTrackHypothesis*> dKinFitChargedHypoMap;
		unordered_map<const DKinFitParticle*, DNeutralParticleHypothesis*> dKinFitNeutralHypoMap;
		unordered_map<const DKinFitParticle*, DBeamPhoton*> dKinFitBeamPhotonMap;

		//RESOURCE POOLS
		DResourcePool<DEventRFBunch> dResourcePool_EventRFBunch;
		DResourcePool<DParticleCombo> dResourcePool_ParticleCombo;
		DResourcePool<DParticleComboStep> dResourcePool_ParticleComboStep;
		DResourcePool<DKinematicData> dResourcePool_KinematicData;
};

}

#endif // DParticleComboCreator_h
