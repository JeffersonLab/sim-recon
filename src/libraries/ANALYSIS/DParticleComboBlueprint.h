#ifndef _DParticleComboBlueprint_
#define _DParticleComboBlueprint_

#include <deque>

#include "JANA/JObject.h"
#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DParticleComboBlueprintStep.h"

using namespace std;
using namespace jana;

class DParticleComboBlueprint : public JObject
{
	public:

		JOBJECT_PUBLIC(DParticleComboBlueprint);

		DParticleComboBlueprint(void) : dReaction(NULL) {}

		const DParticleComboBlueprintStep* Get_ParticleComboBlueprintStep(size_t locStepIndex) const;
		const DParticleComboBlueprintStep* Pop_ParticleComboBlueprintStep(void);
		inline void Add_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep){dParticleComboBlueprintSteps.push_back(locParticleComboBlueprintStep);}
		inline size_t Get_NumParticleComboBlueprintSteps(void) const{return dParticleComboBlueprintSteps.size();}
		void Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locStepIndex);

		inline const DReaction* Get_Reaction(void) const{return dReaction;}
		inline void Set_Reaction(const DReaction* locReaction){dReaction = locReaction;}

	private:
		const DReaction* dReaction;
		deque<const DParticleComboBlueprintStep*> dParticleComboBlueprintSteps; //must be in order you want to evaluate them
};

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
	const DParticleComboBlueprintStep* locParticleComboBlueprintStep = dParticleComboBlueprintSteps.back();
	dParticleComboBlueprintSteps.pop_back();
	return locParticleComboBlueprintStep;
}

#endif // _DParticleComboBlueprint_

