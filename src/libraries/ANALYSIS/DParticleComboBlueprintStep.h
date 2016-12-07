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
		bool operator<(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const;
		inline bool operator!=(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const{return (!((*this) == locParticleComboBlueprintStep));}
		void Reset(void);

		inline const DReactionStep* Get_ReactionStep(void) const{return dReactionStep;}
		inline void Set_ReactionStep(const DReactionStep* locReactionStep){dReactionStep = locReactionStep;}

		void Add_FinalParticle_SourceObject(const JObject* locObject, int locDecayStepIndex);
		const JObject* Pop_FinalParticle_SourceObject(void);
		inline size_t Get_NumFinalParticleSourceObjects(void) const{return dFinalParticleSourceObjects.size();}
		const JObject* Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const;

		void Get_DetectedNeutralShowerSourceObjects(deque<pair<const DNeutralShower*, Particle_t> >& locNeutralShowers) const;
		void Get_DetectedChargedTrackSourceObjects(deque<pair<const DChargedTrack*, Particle_t> >& locChargedTracks) const;

		inline Particle_t Get_InitialParticleID(void) const{return ((dReactionStep != NULL) ? dReactionStep->Get_InitialParticleID() : Unknown);}
		inline Particle_t Get_TargetParticleID(void) const{return ((dReactionStep != NULL) ? dReactionStep->Get_TargetParticleID() : Unknown);}

		inline int Get_InitialParticleDecayFromStepIndex(void) const{return int(dInitialParticleDecayFromStepIndex);}
		inline void Set_InitialParticleDecayFromStepIndex(int locInitialParticleDecayFromStepIndex){dInitialParticleDecayFromStepIndex = (signed char)locInitialParticleDecayFromStepIndex;}

		inline Particle_t Get_FinalParticleID(size_t locFinalParticleIndex) const{return ((dReactionStep != NULL) ? dReactionStep->Get_FinalParticleID(locFinalParticleIndex) : Unknown);}
		void Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const;

		inline size_t Get_NumDecayStepIndices(void) const{return dDecayStepIndices.size();}
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
		inline bool Is_InitialParticleDetected(void) const{return (Get_InitialParticleID() == Gamma);}
		inline bool Is_InitialParticleDecaying(void) const{return (Get_InitialParticleID() != Gamma);}
		inline bool Is_InitialParticleMissing(void) const{return false;} //currently not supported!
		inline bool Is_InitialParticleCharged(void) const{return (ParticleCharge(Get_InitialParticleID()) != 0);}
		inline bool Is_InitialParticleNeutral(void) const{return (ParticleCharge(Get_InitialParticleID()) == 0);}

		//target particle
		inline bool Is_TargetPresent(void) const{return (Get_TargetParticleID() != Unknown);}
		inline bool Is_TargetParticleCharged(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) != 0) : false);}
		inline bool Is_TargetParticleNeutral(void) const{return ((Get_TargetParticleID() != Unknown) ? (ParticleCharge(Get_TargetParticleID()) == 0) : false);}

	private:
		const DReactionStep* dReactionStep;

		deque<const JObject*> dFinalParticleSourceObjects; //NULL if decaying or missing
		signed char dInitialParticleDecayFromStepIndex; //-1 if photon, else index points to step index it is produced at
		deque<signed char> dDecayStepIndices; //one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
};

struct DParticleComboBlueprintStep_Comparator
{
	bool operator() (DParticleComboBlueprintStep* const& lhs, DParticleComboBlueprintStep* const& rhs) const
	{
		return (*lhs) < (*rhs);
	}
};

inline void DParticleComboBlueprintStep::Reset(void)
{
	dReactionStep = NULL;
	dInitialParticleDecayFromStepIndex = -1;
	dFinalParticleSourceObjects.clear();
	dDecayStepIndices.clear();
}

inline void DParticleComboBlueprintStep::Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const
{
	if(dReactionStep != NULL)
		dReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
}

inline void DParticleComboBlueprintStep::Add_FinalParticle_SourceObject(const JObject* locObject, int locDecayStepIndex)
{
	dFinalParticleSourceObjects.push_back(locObject);
	dDecayStepIndices.push_back((signed char)locDecayStepIndex);
}

inline const JObject* DParticleComboBlueprintStep::Pop_FinalParticle_SourceObject(void)
{
	if(dFinalParticleSourceObjects.empty())
		return NULL;
	const JObject* locObject = dFinalParticleSourceObjects.back();
	dFinalParticleSourceObjects.pop_back();
	dDecayStepIndices.pop_back();
	return locObject;
}

inline const JObject* DParticleComboBlueprintStep::Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticleSourceObjects.size())
		return NULL;
	return dFinalParticleSourceObjects[locFinalParticleIndex];
}

inline int DParticleComboBlueprintStep::Get_DecayStepIndex(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dDecayStepIndices.size())
		return -1;
	return int(dDecayStepIndices[locFinalParticleIndex]);
}

inline int DParticleComboBlueprintStep::Get_MissingParticleIndex(void) const //-1 for no missing particles, else final state particle at this index is missing
{
	for(size_t loc_i = 0; loc_i < dDecayStepIndices.size(); ++loc_i)
	{
		if(dDecayStepIndices[loc_i] == ((signed char)-1))
			return loc_i;
	}
	return -1;
}

inline bool DParticleComboBlueprintStep::Is_FinalParticleCharged(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticleSourceObjects())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) != 0);
}

inline bool DParticleComboBlueprintStep::Is_FinalParticleNeutral(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticleSourceObjects())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) == 0);
}

inline void DParticleComboBlueprintStep::Get_DetectedNeutralShowerSourceObjects(deque<pair<const DNeutralShower*, Particle_t> >& locNeutralShowers) const
{
	locNeutralShowers.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleSourceObjects.size(); ++loc_i)
	{
		if(dFinalParticleSourceObjects[loc_i] == NULL)
			continue;
		Particle_t locPID = Get_FinalParticleID(loc_i);
		if(ParticleCharge(locPID) != 0)
			continue;
		pair<const DNeutralShower*, Particle_t> locPIDPair(dynamic_cast<const DNeutralShower*>(dFinalParticleSourceObjects[loc_i]), locPID);
		locNeutralShowers.push_back(locPIDPair);
	}
}

inline void DParticleComboBlueprintStep::Get_DetectedChargedTrackSourceObjects(deque<pair<const DChargedTrack*, Particle_t> >& locChargedTracks) const
{
	locChargedTracks.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticleSourceObjects.size(); ++loc_i)
	{
		if(dFinalParticleSourceObjects[loc_i] == NULL)
			continue;
		Particle_t locPID = Get_FinalParticleID(loc_i);
		if(ParticleCharge(locPID) == 0)
			continue;
		pair<const DChargedTrack*, Particle_t> locPIDPair(dynamic_cast<const DChargedTrack*>(dFinalParticleSourceObjects[loc_i]), locPID);
		locChargedTracks.push_back(locPIDPair);
	}
}

#endif // _DParticleComboBlueprintStep_

