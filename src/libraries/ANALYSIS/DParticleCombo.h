#ifndef _DParticleCombo_
#define _DParticleCombo_

#include <deque>
#include <set>

#include "JANA/JObject.h"
#include "particleType.h"
#include "PID/DKinematicData.h"
#include "PID/DEventRFBunch.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DKinFitResults.h"

using namespace std;
using namespace jana;

class DParticleCombo : public JObject
{
	public:
		JOBJECT_PUBLIC(DParticleCombo);

		// SET STEPS
		inline void Add_ParticleComboStep(const DParticleComboStep* locParticleComboStep){dParticleComboSteps.push_back(locParticleComboStep);}
		void Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex);

		// SET OBJECT DATA:
		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}
		inline void Set_KinFitResults(const DKinFitResults* locKinFitResults){dKinFitResults = locKinFitResults;}
		inline void Set_EventRFBunch(const DEventRFBunch* locEventRFBunch){dEventRFBunch = locEventRFBunch;}

		// GET OBJECT DATA:
		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		inline const DKinFitResults* Get_KinFitResults(void) const{return dKinFitResults;}
		inline const DEventRFBunch* Get_EventRFBunch(void) const{return dEventRFBunch;}

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

		// GET FINAL PARTICLE SOURCE OBJECTS - BY TRAIT:
		void Get_DetectedFinalChargedParticles_SourceObjects(deque<const DChargedTrack*>& locSourceChargedTracks) const;
		void Get_DetectedFinalNeutralParticles_SourceObjects(deque<const DNeutralShower*>& locSourceNeutralShowers) const;

		// OTHER:
		set<pair<const JObject*, Particle_t> > Get_DecayingParticleSourceObjects(size_t locStepIndex) const;
		bool Get_ApplyKinFitMassConstraintOnInitialParticleFlag(size_t locStepIndex) const;
		bool Check_AreMeasuredParticlesIdentical(const DParticleCombo* locParticleCombo) const;

	private:
		// PRIVATE METHODS:
		void Get_DecayChainParticles_Measured_Recursive(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const;

		const DReaction* dReaction;
		const DKinFitResults* dKinFitResults;
		const DEventRFBunch* dEventRFBunch;
		deque<const DParticleComboStep*> dParticleComboSteps;
};

inline const DParticleComboStep* DParticleCombo::Get_ParticleComboStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboSteps.size())
		return NULL;
	return dParticleComboSteps[locStepIndex];
}

inline void DParticleCombo::Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex)
{
	if(locStepIndex >= Get_NumParticleComboSteps())
		return;
	dParticleComboSteps[locStepIndex] = locParticleComboStep;
}

inline bool DParticleCombo::Get_ApplyKinFitMassConstraintOnInitialParticleFlag(size_t locStepIndex) const
{
	if(dReaction == NULL)
		return false;
	return dReaction->Get_ReactionStep(locStepIndex)->Get_ApplyKinFitMassConstraintOnInitialParticleFlag();
}

inline void DParticleCombo::Get_ParticleComboSteps(Particle_t locInitialPID, deque<const DParticleComboStep*>& locParticleComboStepDeque) const
{
	locParticleComboStepDeque.clear();
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		if(dParticleComboSteps[loc_i]->Get_InitialParticleID() == locInitialPID)
			locParticleComboStepDeque.push_back(dParticleComboSteps[loc_i]);
	}
}

inline void DParticleCombo::Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline const DKinematicData* DParticleCombo::Get_MissingParticle(void) const
{
	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locKinematicData = dParticleComboSteps[loc_i]->Get_MissingParticle();
		if(locKinematicData != NULL)
			return locKinematicData;
	}
	return NULL;
}

inline void DParticleCombo::Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_FinalParticles(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_FinalParticles_Measured(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

inline void DParticleCombo::Get_DecayChainParticles_Measured(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const
{
	locMeasuredParticles.clear();
	Get_DecayChainParticles_Measured_Recursive(locStepIndex, locMeasuredParticles);
}

#endif // _DParticleCombo_
