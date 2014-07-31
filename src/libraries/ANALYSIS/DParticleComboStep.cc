#include "DParticleComboStep.h"

bool DParticleComboStep::operator==(const DParticleComboStep& locParticleComboStep) const
{
	if(dParticleComboBlueprintStep != locParticleComboStep.dParticleComboBlueprintStep)
		return false;

	if(dInitialParticle != locParticleComboStep.dInitialParticle)
		return false;
	if(dInitialParticle_Measured != locParticleComboStep.dInitialParticle_Measured)
		return false;
	if(dTargetParticle != locParticleComboStep.dTargetParticle)
		return false;

	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(dFinalParticles[loc_i] != locParticleComboStep.dFinalParticles[loc_i])
			return false;
	}

	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
	{
		if(dFinalParticles_Measured[loc_i] != locParticleComboStep.dFinalParticles_Measured[loc_i])
			return false;
	}
	return true;
}

string DParticleComboStep::Get_FinalParticlesROOTName(void) const
{
	deque<Particle_t> locFinalParticleIDs;
	if(dParticleComboBlueprintStep != NULL)
		dParticleComboBlueprintStep->Get_FinalParticleIDs(locFinalParticleIDs);

	string locStepROOTName;
	for(size_t loc_i = 0; loc_i < locFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == Get_MissingParticleIndex())
			continue;
		locStepROOTName += ParticleName_ROOT(locFinalParticleIDs[loc_i]);
	}
	if(Get_MissingParticleIndex() >= 0)
		locStepROOTName += string("(") + ParticleName_ROOT(locFinalParticleIDs[Get_MissingParticleIndex()]) + string(")");
	return locStepROOTName;
}

void DParticleComboStep::Get_FinalParticlesROOTName(deque<string>& locParticleNames) const
{
	deque<Particle_t> locFinalParticleIDs;
	if(dParticleComboBlueprintStep != NULL)
		dParticleComboBlueprintStep->Get_FinalParticleIDs(locFinalParticleIDs);

	locParticleNames.clear();
	for(size_t loc_i = 0; loc_i < locFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == Get_MissingParticleIndex())
			continue;
		locParticleNames.push_back(ParticleName_ROOT(locFinalParticleIDs[loc_i]));
	}
	if(Get_MissingParticleIndex() >= 0)
		locParticleNames.push_back(string("(") + ParticleName_ROOT(locFinalParticleIDs[Get_MissingParticleIndex()]) + string(")"));
}

string DParticleComboStep::Get_FinalDetectedParticlesROOTName(void) const
{
	deque<Particle_t> locFinalParticleIDs;
	if(dParticleComboBlueprintStep != NULL)
		dParticleComboBlueprintStep->Get_FinalParticleIDs(locFinalParticleIDs);

	string locStepROOTName;
	for(size_t loc_i = 0; loc_i < locFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == Get_MissingParticleIndex())
			continue;
		locStepROOTName += ParticleName_ROOT(locFinalParticleIDs[loc_i]);
	}
	return locStepROOTName;
}

string DParticleComboStep::Get_StepName(void) const
{
	string locStepName = ParticleType(Get_InitialParticleID());
	if(Get_TargetParticleID() != Unknown)
		locStepName += string("_") + ParticleType(Get_TargetParticleID());
	locStepName += "_->";

	deque<Particle_t> locFinalParticleIDs;
	if(dParticleComboBlueprintStep != NULL)
		dParticleComboBlueprintStep->Get_FinalParticleIDs(locFinalParticleIDs);

	for(size_t loc_i = 0; loc_i < locFinalParticleIDs.size(); ++loc_i)
	{
		if(int(loc_i) == Get_MissingParticleIndex())
			continue;
		locStepName += string("_") + ParticleType(locFinalParticleIDs[loc_i]);
	}
	if(Get_MissingParticleIndex() >= 0)
		locStepName += string("_(") + ParticleType(locFinalParticleIDs[Get_MissingParticleIndex()]) + string(")");
	return locStepName;
}


