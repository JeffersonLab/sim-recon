#include "ANALYSIS/DParticleComboBlueprintStep.h"

bool DParticleComboBlueprintStep::operator<(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const
{
	if(dReactionStep < locParticleComboBlueprintStep.dReactionStep)
		return true;
	else if(dReactionStep > locParticleComboBlueprintStep.dReactionStep)
		return false;

	if(dInitialParticleDecayFromStepIndex < locParticleComboBlueprintStep.dInitialParticleDecayFromStepIndex)
		return true;
	else if(dInitialParticleDecayFromStepIndex > locParticleComboBlueprintStep.dInitialParticleDecayFromStepIndex)
		return false;

	if(dFinalParticleSourceObjects.size() < locParticleComboBlueprintStep.dFinalParticleSourceObjects.size())
		return true;
	if(dFinalParticleSourceObjects.size() > locParticleComboBlueprintStep.dFinalParticleSourceObjects.size())
		return false;

	for(size_t loc_i = 0; loc_i < dFinalParticleSourceObjects.size(); ++loc_i)
	{
		if(dFinalParticleSourceObjects[loc_i] < locParticleComboBlueprintStep.dFinalParticleSourceObjects[loc_i])
			return true;
		if(dFinalParticleSourceObjects[loc_i] > locParticleComboBlueprintStep.dFinalParticleSourceObjects[loc_i])
			return false;
	}

	if(dDecayStepIndices.size() < locParticleComboBlueprintStep.dDecayStepIndices.size())
		return true;
	if(dDecayStepIndices.size() > locParticleComboBlueprintStep.dDecayStepIndices.size())
		return false;

	for(size_t loc_i = 0; loc_i < dDecayStepIndices.size(); ++loc_i)
	{
		if(dDecayStepIndices[loc_i] < locParticleComboBlueprintStep.dDecayStepIndices[loc_i])
			return true;
		else if(dDecayStepIndices[loc_i] > locParticleComboBlueprintStep.dDecayStepIndices[loc_i])
			return false;
	}

	return false; //equivalent!
}

bool DParticleComboBlueprintStep::operator==(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const
{
	if(dReactionStep != locParticleComboBlueprintStep.dReactionStep)
		return false;

	if(dInitialParticleDecayFromStepIndex != locParticleComboBlueprintStep.dInitialParticleDecayFromStepIndex)
		return false;

	if(dFinalParticleSourceObjects.size() != locParticleComboBlueprintStep.dFinalParticleSourceObjects.size())
		return false;

	for(size_t loc_i = 0; loc_i < dFinalParticleSourceObjects.size(); ++loc_i)
	{
		if(dFinalParticleSourceObjects[loc_i] != locParticleComboBlueprintStep.dFinalParticleSourceObjects[loc_i])
			return false;
	}

	if(dDecayStepIndices.size() != locParticleComboBlueprintStep.dDecayStepIndices.size())
		return false;

	for(size_t loc_i = 0; loc_i < dDecayStepIndices.size(); ++loc_i)
	{
		if(dDecayStepIndices[loc_i] != locParticleComboBlueprintStep.dDecayStepIndices[loc_i])
			return false;
	}

	return true;
}

