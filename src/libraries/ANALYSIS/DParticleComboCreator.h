#ifndef DParticleComboCreator_h
#define DParticleComboCreator_h

#include <unordered_map>
#include <map>

#include "JANA/JEventLoop.h"
#include "PID/DEventRFBunch.h"
#include "PID/DChargedTrackHypothesis_factory.h"
#include "PID/DNeutralParticleHypothesis_factory.h"
#include "PID/DBeamPhoton_factory.h"
#include "PID/DParticleID.h"

#include "ANALYSIS/DKinFitUtils_GlueX.h"
#include "ANALYSIS/DParticleCombo.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DSourceComboTimeHandler.h"
#include "ANALYSIS/DSourceComboVertexer.h"

using namespace std;
using namespace jana;

class DAnalysisUtilities;

namespace DAnalysis
{

class DSourceComboer;

class DParticleComboCreator
{
	public:
		DParticleComboCreator(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboTimeHandler* locSourceComboTimeHandler, const DSourceComboVertexer* locSourceComboVertexer);

		const DParticleCombo* Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift, DKinFitType locKinFitType);
		const DParticleCombo* Create_KinFitCombo_NewCombo(const DParticleCombo* locOrigCombo, const DReaction* locReaction, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain);

		const DParticleCombo* Build_ThrownCombo(JEventLoop* locEventLoop);
		const DParticleCombo* Build_ThrownCombo(JEventLoop* locEventLoop, const DReaction* locThrownReaction, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps);

		void Reset(void);
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}

	private:

		int dDebugLevel = 0;

		bool Get_CreateNeutralErrorMatrixFlag_Combo(const DReactionVertexInfo* locReactionVertexInfo, DKinFitType locKinFitType);

		//DECAYING PARTICLES, POST-KINFIT
		void Set_DecayingParticles(const DReaction* locReaction, const DParticleCombo* locNewParticleCombo, const DParticleCombo* locOldParticleCombo, size_t locStepIndex, DParticleComboStep* locNewParticleComboStep, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		DKinFitParticle* Get_DecayingParticle(const DReaction* locReaction, const DParticleCombo* locOldParticleCombo, size_t locComboStepIndex, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults);
		bool Search_ForParticleInDecay(const DKinFitChain* locKinFitChain, size_t locStepToSearch, DKinFitParticle* locParticleToFind);

		//SPACETIME VERTEX, POST-KINFIT
		void Set_SpacetimeVertex(const DReaction* locReaction, const DParticleCombo* locNewParticleCombo, DParticleComboStep* locNewParticleComboStep, size_t locStepIndex, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain) const;

		//CREATE PARTICLES
		const DChargedTrackHypothesis* Create_ChargedHypo(const DChargedTrack* locChargedTrack, Particle_t locPID, double locPropagatedRFTime, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryFullCombo, const DKinematicData* locBeamParticle);
		const DBeamPhoton* Create_BeamPhoton_KinFit(const DBeamPhoton* locBeamPhoton, const DKinFitParticle* locKinFitParticle, const DLorentzVector& locSpacetimeVertex);
		const DChargedTrackHypothesis* Create_ChargedHypo_KinFit(const DChargedTrack* locChargedTrack, Particle_t locPID, const DKinFitParticle* locKinFitParticle);
		const DNeutralParticleHypothesis* Create_NeutralHypo_KinFit(const DNeutralParticleHypothesis* locOrigHypo, DKinFitParticle* locKinFitParticle);
		DKinematicData* Build_KinematicData(DKinFitParticle* locKinFitParticle, DKinFitType locKinFitType, DVector3 locPreKinFitVertex);

		TMatrixFSym dVertexCovMatrix;
		unordered_map<const DReactionVertexInfo*, bool> dDanglingNeutralsFlagMap;

		//UTILITIES
		const DSourceComboer* dSourceComboer = nullptr;
		const DSourceComboTimeHandler* dSourceComboTimeHandler = nullptr;
		const DSourceComboVertexer* dSourceComboVertexer = nullptr;
		const DParticleID* dParticleID = nullptr;
		const DAnalysisUtilities* dAnalysisUtilities = nullptr;
		DKinFitUtils_GlueX* dKinFitUtils = nullptr;

		//FACTORIES
		DNeutralParticleHypothesis_factory* dNeutralParticleHypothesisFactory;
		DChargedTrackHypothesis_factory* dChargedTrackHypothesisFactory;
		DBeamPhoton_factory* dBeamPhotonfactory;

		//CREATED OBJECT MAPS
		map<tuple<DReactionStep, const DSourceCombo*, bool, bool, const DSourceCombo*, const DKinematicData*>, const DParticleComboStep*> dComboStepMap; //kindata is beam (null if not in step): for vertex
		unordered_map<int, const DEventRFBunch*> dRFBunchMap;
		map<tuple<const DChargedTrack*, Particle_t, int, bool, const DSourceCombo*, const DSourceCombo*, const DKinematicData*>, const DChargedTrackHypothesis*> dChargedHypoMap;
		map<tuple<const DNeutralShower*, Particle_t, int, bool, bool, const DSourceCombo*, const DSourceCombo*, const DKinematicData*>, const DNeutralParticleHypothesis*> dNeutralHypoMap;
		map<tuple<const DReactionVertexInfo*, const DSourceCombo*, const DKinematicData*, int, bool>, DParticleCombo*> dComboMap;
		unordered_map<const DKinFitParticle*, DChargedTrackHypothesis*> dKinFitChargedHypoMap;
		unordered_map<const DKinFitParticle*, DNeutralParticleHypothesis*> dKinFitNeutralHypoMap;
		unordered_map<const DKinFitParticle*, DBeamPhoton*> dKinFitBeamPhotonMap;

		//RESOURCE POOLS
		DResourcePool<DEventRFBunch> dResourcePool_EventRFBunch;
		DResourcePool<DParticleCombo> dResourcePool_ParticleCombo;
		DResourcePool<DParticleComboStep> dResourcePool_ParticleComboStep;
		DResourcePool<DKinematicData> dResourcePool_KinematicData;
		vector<DKinematicData*> dCreated_KinematicData;
		vector<DParticleCombo*> dCreated_ParticleCombo;
		vector<DParticleComboStep*> dCreated_ParticleComboStep;
		vector<DChargedTrackHypothesis*> dCreated_ChargedHypo;
		vector<DNeutralParticleHypothesis*> dCreated_NeutralHypo;
		vector<DBeamPhoton*> dCreated_BeamPhoton;

		DParticleCombo* Get_ParticleComboResource(void)
		{
			auto locResource = dResourcePool_ParticleCombo.Get_Resource();
			locResource->Reset();
			dCreated_ParticleCombo.push_back(locResource);
			return locResource;
		}

		DParticleComboStep* Get_ParticleComboStepResource(void)
		{
			auto locResource = dResourcePool_ParticleComboStep.Get_Resource();
			locResource->Reset();
			dCreated_ParticleComboStep.push_back(locResource);
			return locResource;
		}
};

}

#endif // DParticleComboCreator_h
