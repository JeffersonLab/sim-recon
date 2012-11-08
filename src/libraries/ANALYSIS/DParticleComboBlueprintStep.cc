#include "ANALYSIS/DParticleComboBlueprintStep.h"

void DParticleComboBlueprintStep::Reset(void)
{
	dReactionStep = NULL;
	dInitialParticleDecayFromStepIndex = -1;
	dFinalParticleSourceObjects.clear();
	dDecayStepIndices.clear();
}

bool DParticleComboBlueprintStep::operator==(const DParticleComboBlueprintStep& locParticleComboBlueprintStep) const
{
	if(dReactionStep != locParticleComboBlueprintStep.dReactionStep)
		return false;

	if(dInitialParticleDecayFromStepIndex != locParticleComboBlueprintStep.dInitialParticleDecayFromStepIndex)
		return false;

	for(size_t loc_i = 0; loc_i < dDecayStepIndices.size(); ++loc_i)
	{
		if(dDecayStepIndices[loc_i] != locParticleComboBlueprintStep.dDecayStepIndices[loc_i])
			return false;
	}
	for(size_t loc_i = 0; loc_i < dFinalParticleSourceObjects.size(); ++loc_i)
	{
		if(dFinalParticleSourceObjects[loc_i] != locParticleComboBlueprintStep.dFinalParticleSourceObjects[loc_i])
			return false;
	}

	return true;
}

void DParticleComboBlueprintStep::Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const
{
	if(dReactionStep != NULL)
		dReactionStep->Get_FinalParticleIDs(locFinalParticleIDs);
}

void DParticleComboBlueprintStep::Add_FinalParticle_SourceObject(const JObject* locObject, int locDecayStepIndex)
{
	dFinalParticleSourceObjects.push_back(locObject);
	dDecayStepIndices.push_back(locDecayStepIndex);
}

const JObject* DParticleComboBlueprintStep::Pop_FinalParticle_SourceObject(void)
{
	if(dFinalParticleSourceObjects.empty())
		return NULL;
	const JObject* locObject = dFinalParticleSourceObjects.back();
	dFinalParticleSourceObjects.pop_back();
	dDecayStepIndices.pop_back();
	return locObject;
}

const JObject* DParticleComboBlueprintStep::Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticleSourceObjects.size())
		return NULL;
	return dFinalParticleSourceObjects[locFinalParticleIndex];
}

int DParticleComboBlueprintStep::Get_DecayStepIndex(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dDecayStepIndices.size())
		return -1;
	return dDecayStepIndices[locFinalParticleIndex];
}

int DParticleComboBlueprintStep::Get_MissingParticleIndex(void) const //-1 for no missing particles, else final state particle at this index is missing
{
	for(size_t loc_i = 0; loc_i < dDecayStepIndices.size(); ++loc_i)
	{
		if(dDecayStepIndices[loc_i] == -1)
			return loc_i;
	}
	return -1;
}

bool DParticleComboBlueprintStep::Is_FinalParticleCharged(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticleSourceObjects())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) != 0);
}

bool DParticleComboBlueprintStep::Is_FinalParticleNeutral(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticleSourceObjects())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) == 0);
}

