#include "ANALYSIS/DParticleComboBlueprint.h"

const DParticleComboBlueprintStep* DParticleComboBlueprint::Get_ParticleComboBlueprintStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboBlueprintSteps.size())
		return NULL;
	return dParticleComboBlueprintSteps[locStepIndex];
}

void DParticleComboBlueprint::Set_ParticleComboBlueprintStep(const DParticleComboBlueprintStep* locParticleComboBlueprintStep, size_t locStepIndex)
{
	if(locStepIndex < dParticleComboBlueprintSteps.size())
		dParticleComboBlueprintSteps[locStepIndex] = locParticleComboBlueprintStep;
}

const DParticleComboBlueprintStep* DParticleComboBlueprint::Pop_ParticleComboBlueprintStep(void)
{
	if(dParticleComboBlueprintSteps.empty())
		return NULL;
	const DParticleComboBlueprintStep* locParticleComboBlueprintStep = dParticleComboBlueprintSteps.back();
	dParticleComboBlueprintSteps.pop_back();
	return locParticleComboBlueprintStep;
}

