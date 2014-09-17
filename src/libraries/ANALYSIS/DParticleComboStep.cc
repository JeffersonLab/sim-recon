#include "DParticleComboStep.h"

DParticleComboStep::DParticleComboStep(void)
{
	Reset();
}

void DParticleComboStep::Reset(void)
{
	dParticleComboBlueprintStep = NULL;
	dInitialParticle = NULL;
	dInitialParticle_Measured = NULL;
	dTargetParticle = NULL;

	dFinalParticles.clear();
	dFinalParticles_Measured.clear();
}

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

void DParticleComboStep::Get_FinalParticleIDs(deque<Particle_t>& locFinalParticleIDs) const
{
	if(dParticleComboBlueprintStep != NULL)
		dParticleComboBlueprintStep->Get_FinalParticleIDs(locFinalParticleIDs);
}

string DParticleComboStep::Get_InitialParticlesROOTName(void) const
{
	string locStepROOTName = ParticleName_ROOT(Get_InitialParticleID());
	if(Get_TargetParticleID() != Unknown)
		locStepROOTName += ParticleName_ROOT(Get_TargetParticleID());
	return locStepROOTName;
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

string DParticleComboStep::Get_StepROOTName(void) const
{
	string locStepROOTName = Get_InitialParticlesROOTName();
	locStepROOTName += "#rightarrow";
	locStepROOTName += Get_FinalParticlesROOTName();
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

const JObject* DParticleComboStep::Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const
{
	if(dParticleComboBlueprintStep == NULL)
		return NULL;
	return dParticleComboBlueprintStep->Get_FinalParticle_SourceObject(locFinalParticleIndex);
}

const DKinematicData* DParticleComboStep::Get_FinalParticle(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticles.size())
		return NULL;
	return dFinalParticles[locFinalParticleIndex];
}

const DKinematicData* DParticleComboStep::Get_FinalParticle_Measured(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticles_Measured.size())
		return NULL;
	return dFinalParticles_Measured[locFinalParticleIndex];
}

void DParticleComboStep::Get_FinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		locParticles.push_back(dFinalParticles[loc_i]);
}

void DParticleComboStep::Get_FinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
		locParticles.push_back(dFinalParticles_Measured[loc_i]);
}

void DParticleComboStep::Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(Get_FinalParticleID(loc_i) == locPID)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

void DParticleComboStep::Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
	{
		if(Get_FinalParticleID(loc_i) == locPID)
			locParticles.push_back(dFinalParticles_Measured[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(Is_FinalParticleDetected(loc_i))
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) == 0)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) != 0)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
	{
		if(Is_FinalParticleDetected(loc_i))
			locParticles.push_back(dFinalParticles_Measured[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) == 0)
			locParticles.push_back(dFinalParticles_Measured[loc_i]);
	}
}

void DParticleComboStep::Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles_Measured.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) != 0)
			locParticles.push_back(dFinalParticles_Measured[loc_i]);
	}
}

const DKinematicData* DParticleComboStep::Get_MissingParticle(void) const
{
	int locMissingParticleIndex = Get_MissingParticleIndex();
	if(locMissingParticleIndex == -1)
		return NULL;
	return dFinalParticles[locMissingParticleIndex];
}

bool DParticleComboStep::Is_FinalParticleCharged(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticles())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) != 0);
}

bool DParticleComboStep::Is_FinalParticleNeutral(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= Get_NumFinalParticles())
		return false;
	return (ParticleCharge(Get_FinalParticleID(locFinalParticleIndex)) == 0);
}

