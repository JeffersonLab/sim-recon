#ifndef _DParticleComboStep_
#define _DParticleComboStep_

#include <deque>
#include <iostream>

#include "JANA/JObject.h"
#include "particleType.h"
#include "PID/DKinematicData.h"
#include "DParticleComboBlueprint.h"

using namespace std;
using namespace jana;

class DParticleComboStep
{
	public:

		// CONSTRUCTOR:
		DParticleComboStep(void);

		// RESET:
		void Reset(void);

		// SET INITIAL AND TARGET PARTICLES:
		inline void Set_InitialParticle(const DKinematicData* locInitialParticle){dInitialParticle = locInitialParticle;}
		inline void Set_InitialParticle_Measured(const DKinematicData* locInitialParticle){dInitialParticle_Measured = locInitialParticle;}
		inline void Set_TargetParticle(const DKinematicData* locTargetParticle){dTargetParticle = locTargetParticle;}

		// SET FINAL PARTICLES:
		inline void Add_FinalParticle(const DKinematicData* locFinalParticle){dFinalParticles.push_back(locFinalParticle);}
		inline void Add_FinalParticle_Measured(const DKinematicData* locFinalParticle){dFinalParticles_Measured.push_back(locFinalParticle);}

		// SET BLUEPRINT:
		inline void Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintStep = locParticleComboBlueprintStep;}

		// GET PARTICLE PIDs:
		inline Particle_t Get_InitialParticleID(void) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_InitialParticleID() : Unknown);}
		inline Particle_t Get_TargetParticleID(void) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_TargetParticleID() : Unknown);}
		inline Particle_t Get_FinalParticleID(size_t locFinalParticleIndex) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_FinalParticleID(locFinalParticleIndex) : Unknown);}
		void Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const;

		// GET INITIAL AND TARGET PARTICLES:
		inline const DKinematicData* Get_InitialParticle(void) const{return dInitialParticle;}
		const DKinematicData* Get_InitialParticle_Measured(void) const{return dInitialParticle_Measured;};
		inline const DKinematicData* Get_TargetParticle(void) const{return dTargetParticle;}

		// GET FINAL PARTICLES:
		inline size_t Get_NumFinalParticles(void) const{return dFinalParticles.size();}
		const DKinematicData* Get_FinalParticle(size_t locFinalParticleIndex) const;
		const DKinematicData* Get_FinalParticle_Measured(size_t locFinalParticleIndex) const;
		const JObject* Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const; //if missing or decaying: source object is null!
		void Get_FinalParticles(deque<const DKinematicData*>& locParticles) const;
		const DKinematicData* Get_MissingParticle(void) const; //returns NULL if none missing!

		// GET FINAL PARTICLES - BY PID:
		void Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const; //get all final particles of the given PID
		void Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const; //get all measured final particles of the given PID

		// GET FINAL PARTICLES - BY TRAIT:
		void Get_FinalParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const;

		// GET TRAITS - INITIAL PARTICLE:
		bool Is_InitialParticleDetected(void) const{return (Get_InitialParticleID() == Gamma);}
		bool Is_InitialParticleDecaying(void) const{return (Get_InitialParticleID() != Gamma);}
		bool Is_InitialParticleMissing(void) const{return false;} //currently not supported!
		bool Is_InitialParticleCharged(void) const{return (ParticleCharge(Get_InitialParticleID()) != 0);}
		bool Is_InitialParticleNeutral(void) const{return (ParticleCharge(Get_InitialParticleID()) == 0);}

		// GET TRAITS - TARGET PARTICLE:
		bool Is_TargetPresent(void) const{return (Get_TargetParticleID() != Unknown);}
		bool Is_TargetParticleCharged(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) != 0) : false);}
		bool Is_TargetParticleNeutral(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) == 0) : false);}

		// GET TRAITS - FINAL PARTICLES:
		inline bool Is_FinalParticleDetected(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) == -2);}
		inline bool Is_FinalParticleDecaying(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) >= 0);}
		inline bool Is_FinalParticleMissing(size_t locFinalParticleIndex) const{return (Get_DecayStepIndex(locFinalParticleIndex) == -1);}
		bool Is_FinalParticleCharged(size_t locFinalParticleIndex) const;
		bool Is_FinalParticleNeutral(size_t locFinalParticleIndex) const;

		// GET PARTICLE NAME STRINGS:
		string Get_StepName(void) const;
		string Get_InitialParticlesROOTName(void) const;
		string Get_FinalParticlesROOTName(void) const;
		void Get_FinalParticlesROOTName(deque<string>& locParticleNames) const;
		string Get_FinalDetectedParticlesROOTName(void) const;
		string Get_StepROOTName(void) const;

		// GET BLUEPRINT
		inline const DParticleComboBlueprintStep* Get_ParticleComboBlueprintStep(void) const{return dParticleComboBlueprintStep;}

		// GET CONTROL VARIABLES:
		//DecayStepIndex: one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
		inline int Get_DecayStepIndex(size_t locFinalParticleIndex) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_DecayStepIndex(locFinalParticleIndex) : -3);}
		//MissingParticleIndex: -1 for no missing particles, else final state particle at this index is missing
		int Get_MissingParticleIndex(void) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_MissingParticleIndex() : -2);}
		//InitialParticleDecayFromStepIndex: -1 if photon, else index points to step index it is produced at
		inline int Get_InitialParticleDecayFromStepIndex(void) const{return ((dParticleComboBlueprintStep != NULL) ? dParticleComboBlueprintStep->Get_InitialParticleDecayFromStepIndex() : -2);}

		// COMPARISON OPERATORS:
		bool operator==(const DParticleComboStep& locParticleComboStep) const;
		bool operator!=(const DParticleComboStep& locParticleComboStep) const{return (!((*this) == locParticleComboStep));}

	private:
		// BLUEPRINT:
		const DParticleComboBlueprintStep* dParticleComboBlueprintStep; //contains PIDs, source objects

		// INITIAL PARTICLES:
		const DKinematicData* dInitialParticle; //if is null: decaying or beam particle not yet set!
		const DKinematicData* dInitialParticle_Measured; //if is null: decaying or beam particle not yet set!
		const DKinematicData* dTargetParticle; //NULL for no target

		// FINAL PARTICLES:
		deque<const DKinematicData*> dFinalParticles; //if particle is null: missing or decaying! //these are DChargedTrackHypothesis or DNeutralParticleHypothesis objects if detected
		deque<const DKinematicData*> dFinalParticles_Measured; //if particle is null: missing or decaying! //these are DChargedTrackHypothesis or DNeutralParticleHypothesis objects if detected
};

#endif // _DParticleComboStep_

