#ifndef _DParticleComboBlueprintStep_
#define _DParticleComboBlueprintStep_

#include <deque>

#include "JANA/JObject.h"
#include "particleType.h"
#include "PID/DKinematicData.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralShower.h"
#include "ANALYSIS/DReactionStep.h"

using namespace std;
using namespace jana;

class DParticleComboBlueprintStep
{
	public:

		bool operator==(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const;
		bool operator!=(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const{return (!((*this) == locParticleComboBlueprintStep));}
		void Reset(void);

		inline const DReactionStep* Get_ReactionStep(void) const{return dReactionStep;}
		inline void Set_ReactionStep(const DReactionStep* locReactionStep){dReactionStep = locReactionStep;}

		void Add_FinalParticle_SourceObject(const JObject* locObject, int locDecayStepIndex);
		const JObject* Pop_FinalParticle_SourceObject(void);
		inline size_t Get_NumFinalParticleSourceObjects(void) const{return dFinalParticleSourceObjects.size();}
		const JObject* Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const;

		inline Particle_t Get_InitialParticleID(void) const{return ((dReactionStep != NULL) ? dReactionStep->Get_InitialParticleID() : Unknown);}
		inline Particle_t Get_TargetParticleID(void) const{return ((dReactionStep != NULL) ? dReactionStep->Get_TargetParticleID() : Unknown);}

		inline int Get_InitialParticleDecayFromStepIndex(void) const{return dInitialParticleDecayFromStepIndex;}
		inline void Set_InitialParticleDecayFromStepIndex(int locInitialParticleDecayFromStepIndex){dInitialParticleDecayFromStepIndex = locInitialParticleDecayFromStepIndex;}

		inline Particle_t Get_FinalParticleID(size_t locFinalParticleIndex) const{return ((dReactionStep != NULL) ? dReactionStep->Get_FinalParticleID(locFinalParticleIndex) : Unknown);}
		void Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const;

		size_t Get_NumDecayStepIndices(void) const{return dDecayStepIndices.size();}
		int Get_DecayStepIndex(size_t locFinalParticleIndex) const;

		int Get_MissingParticleIndex(void) const; //-1 for no missing particles, else final state particle at this index is missing

		//ASKERS

		//final particles
		inline bool Is_FinalParticleDetected(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) == -2);}
		inline bool Is_FinalParticleDecaying(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) >= 0);}
		inline bool Is_FinalParticleMissing(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) == -1);}
		bool Is_FinalParticleCharged(size_t locFinalParticleIndex) const;
		bool Is_FinalParticleNeutral(size_t locFinalParticleIndex) const;

		//initial particle
		bool Is_InitialParticleDetected(void) const{return (Get_InitialParticleID() == Gamma);}
		bool Is_InitialParticleDecaying(void) const{return (Get_InitialParticleID() != Gamma);}
		bool Is_InitialParticleMissing(void) const{return false;} //currently not supported!
		bool Is_InitialParticleCharged(void) const{return (ParticleCharge(Get_InitialParticleID()) != 0);}
		bool Is_InitialParticleNeutral(void) const{return (ParticleCharge(Get_InitialParticleID()) == 0);}

		//target particle
		bool Is_TargetPresent(void) const{return (Get_TargetParticleID() != Unknown);}
		bool Is_TargetParticleCharged(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) != 0) : false);}
		bool Is_TargetParticleNeutral(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) == 0) : false);}

	private:
		const DReactionStep* dReactionStep;

		deque<const JObject*> dFinalParticleSourceObjects; //NULL if decaying or missing
		int dInitialParticleDecayFromStepIndex; //-1 if photon, else index points to step index it is produced at
		deque<int> dDecayStepIndices; //one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
};

#endif // _DParticleComboBlueprintStep_

