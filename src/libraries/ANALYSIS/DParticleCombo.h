#ifndef _DParticleCombo_
#define _DParticleCombo_

#include <deque>

#include "JANA/JObject.h"
#include "particleType.h"
#include "PID/DKinematicData.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DKinFitResults.h"

using namespace std;
using namespace jana;

class DParticleCombo : public JObject
{
	public:
		JOBJECT_PUBLIC(DParticleCombo);

		// SET OBJECT DATA:
		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}
		inline void Set_KinFitResults(const DKinFitResults* locKinFitResults){dKinFitResults = locKinFitResults;}
		inline void Add_ParticleComboStep(const DParticleComboStep* locParticleComboStep){dParticleComboSteps.push_back(locParticleComboStep);}
		void Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex);

		// GET REACTION AND KINFITRESULTS:
		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		inline const DKinFitResults* Get_KinFitResults(void) const{return dKinFitResults;}

		// GET PARTCILE COMBO STEPS:
		inline size_t Get_NumParticleComboSteps(void) const{return dParticleComboSteps.size();}
		const DParticleComboStep* Get_ParticleComboStep(size_t locStepIndex) const;
		void Get_ParticleComboSteps(Particle_t locInitialPID, deque<const DParticleComboStep*>& locParticleComboStepDeque) const;

		// GET FINAL PARTICLES - BY PID:
		void Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const;
		void Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const;
		void Get_FinalParticles(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const;
		void Get_FinalParticles_Measured(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const;

		// GET FINAL PARTICLES - BY TRAIT:
		void Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const;
		const DKinematicData* Get_MissingParticle(void) const;

		// GET FINAL PARTICLES - BY DECAY CHAIN:
		//get all of the measured particles included in the decaychain starting at locStepIndex
		void Get_DecayChainParticles_Measured(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const;

		// GET DECAY CHAIN PARTICLE NAMES:
		string Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, bool locKinFitResultsFlag = false) const;
		string Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag = false) const;

		// OTHER:
		bool Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const;
		bool Will_KinFitBeIdentical(const DParticleCombo* locParticleCombo) const; //the pointers for the steps must be identical for this to be true!!

	private:
		// PRIVATE METHODS:
		string Get_DecayChainFinalParticlesROOTName_Recursive(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag = false) const;

		const DReaction* dReaction;
		const DKinFitResults* dKinFitResults;
		deque<const DParticleComboStep*> dParticleComboSteps;
};

#endif // _DParticleCombo_
