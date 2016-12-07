#ifndef _DParticleComboBlueprint_
#define _DParticleComboBlueprint_

#include <deque>

#include "JANA/JObject.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleComboBlueprintStep.h"

using namespace std;
using namespace jana;

class DParticleComboBlueprint_factory;

class DParticleComboBlueprint : public JObject
{
	friend class DParticleComboBlueprint_factory;

	public:

		JOBJECT_PUBLIC(DParticleComboBlueprint);

		DParticleComboBlueprint(void) : dReaction(NULL) {}
		void Reset(void);

		const DParticleComboBlueprintStep* Get_ParticleComboBlueprintStep(size_t locStepIndex) const;
		const DParticleComboBlueprintStep* Pop_ParticleComboBlueprintStep(void);
		inline void Prepend_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintSteps.push_front(locParticleComboBlueprintStep);}
		inline size_t Get_NumParticleComboBlueprintSteps(void) const{return dParticleComboBlueprintSteps.size();}
		void Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locStepIndex);

		void Get_DetectedNeutralShowerSourceObjects(deque<pair<const DNeutralShower*, Particle_t> >& locNeutralShowers) const;
		void Get_DetectedChargedTrackSourceObjects(deque<pair<const DChargedTrack*, Particle_t> >& locChargedTracks) const;

		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}

	private:
		const DReaction* dReaction;
		deque<const DParticleComboBlueprintStep*> dParticleComboBlueprintSteps; //must be in order you want to evaluate them
};

inline void DParticleComboBlueprint::Reset(void)
{
	dReaction = NULL;
	dParticleComboBlueprintSteps.clear();
}

inline const DParticleComboBlueprintStep* DParticleComboBlueprint::Get_ParticleComboBlueprintStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboBlueprintSteps.size())
		return NULL;
	return dParticleComboBlueprintSteps[locStepIndex];
}

inline void DParticleComboBlueprint::Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locStepIndex)
{
	if(locStepIndex < dParticleComboBlueprintSteps.size())
		dParticleComboBlueprintSteps[locStepIndex] = locParticleComboBlueprintStep;
}

inline const DParticleComboBlueprintStep* DParticleComboBlueprint::Pop_ParticleComboBlueprintStep(void)
{
	if(dParticleComboBlueprintSteps.empty())
		return NULL;
	const DParticleComboBlueprintStep* locParticleComboBlueprintStep = dParticleComboBlueprintSteps.front();
	dParticleComboBlueprintSteps.pop_front();
	return locParticleComboBlueprintStep;
}

inline void DParticleComboBlueprint::Get_DetectedNeutralShowerSourceObjects(deque<pair<const DNeutralShower*, Particle_t> >& locNeutralShowers) const
{
	locNeutralShowers.clear();
	for(size_t loc_i = 0; loc_i < dParticleComboBlueprintSteps.size(); ++loc_i)
	{
		deque<pair<const DNeutralShower*, Particle_t> > locStepNeutralShowers;
		dParticleComboBlueprintSteps[loc_i]->Get_DetectedNeutralShowerSourceObjects(locStepNeutralShowers);
		locNeutralShowers.insert(locNeutralShowers.end(), locStepNeutralShowers.begin(), locStepNeutralShowers.end());
	}
}

inline void DParticleComboBlueprint::Get_DetectedChargedTrackSourceObjects(deque<pair<const DChargedTrack*, Particle_t> >& locChargedTracks) const
{
	locChargedTracks.clear();
	for(size_t loc_i = 0; loc_i < dParticleComboBlueprintSteps.size(); ++loc_i)
	{
		deque<pair<const DChargedTrack*, Particle_t> > locStepChargedTracks;
		dParticleComboBlueprintSteps[loc_i]->Get_DetectedChargedTrackSourceObjects(locStepChargedTracks);
		locChargedTracks.insert(locChargedTracks.end(), locStepChargedTracks.begin(), locStepChargedTracks.end());
	}
}

#endif // _DParticleComboBlueprint_

