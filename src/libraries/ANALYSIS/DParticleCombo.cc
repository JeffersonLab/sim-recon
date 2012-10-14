#include "ANALYSIS/DParticleCombo.h"

bool DParticleCombo::Will_KinFitBeIdentical(const DParticleCombo* locParticleCombo) const
{
	//the pointers for the steps must be identical for this to be true!!
	if(dParticleComboSteps.size() != locParticleCombo->Get_NumParticleComboSteps())
		return false;
	if(dReaction->Get_KinFitType() != locParticleCombo->Get_Reaction()->Get_KinFitType())
		return false;
	for(size_t loc_k = 0; loc_k < locParticleCombo->Get_NumParticleComboSteps(); ++loc_k)
	{
		if(dParticleComboSteps[loc_k] != locParticleCombo->Get_ParticleComboStep(loc_k))
			return false;
	}
	return true;
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, bool locKinFitResultsFlag) const
{
	deque<string> locParticleNames;
	return Get_DecayChainFinalParticlesROOTName(locStepIndex, locParticleNames, locKinFitResultsFlag);
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag) const
{
	locParticleNames.clear();
	return Get_DecayChainFinalParticlesROOTName_Recursive(locStepIndex, locParticleNames, locKinFitResultsFlag);
}

string DParticleCombo::Get_DecayChainFinalParticlesROOTName_Recursive(size_t locStepIndex, deque<string>& locParticleNames, bool locKinFitResultsFlag) const
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances and excluded particles)
	deque<Particle_t> locPIDs;
	dParticleComboSteps[locStepIndex]->Get_FinalParticleIDs(locPIDs);
	string locFullParticleName;
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		if(dParticleComboSteps[locStepIndex]->Is_FinalParticleDecaying(loc_j))
		{
			//expand if: not-kinfitting, not a fixed mass, OR excluded from kinfit
			if((!locKinFitResultsFlag) || (!IsFixedMass(locPIDs[loc_j])) || Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex))
			{
				locFullParticleName += Get_DecayChainFinalParticlesROOTName_Recursive(dParticleComboSteps[locStepIndex]->Get_DecayStepIndex(loc_j), locParticleNames, locKinFitResultsFlag);
				continue;
			}
		}

		locParticleNames.push_back(ParticleName_ROOT(locPIDs[loc_j]));
		locFullParticleName += ParticleName_ROOT(locPIDs[loc_j]);
		continue;
	}
	return locFullParticleName;
}

const DParticleComboStep* DParticleCombo::Get_ParticleComboStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboSteps.size())
		return NULL;
	return dParticleComboSteps[locStepIndex];
}

void DParticleCombo::Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex)
{
	if(locStepIndex >= Get_NumParticleComboSteps())
		return;
	dParticleComboSteps[locStepIndex] = locParticleComboStep;
}

bool DParticleCombo::Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const
{
	if(dReaction == NULL)
		return false;
	return	 dReaction->Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex);
}

void DParticleCombo::Get_ParticleComboSteps(Particle_t locInitialPID, deque<const DParticleComboStep*>& locParticleComboStepDeque) const
{
	locParticleComboStepDeque.clear();
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		if(dParticleComboSteps[loc_i]->Get_InitialParticleID() == locInitialPID)
			locParticleComboStepDeque.push_back(dParticleComboSteps[loc_i]);
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles_Measured(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalChargedParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalNeutralParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_DetectedFinalParticles(locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

const DKinematicData* DParticleCombo::Get_MissingParticle(void) const
{
	const DKinematicData* locKinematicData;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locKinematicData = dParticleComboSteps[loc_i]->Get_MissingParticle();
		if(locKinematicData != NULL)
			return locKinematicData;
	}
	return NULL;
}

void DParticleCombo::Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < dParticleComboSteps.size(); ++loc_i)
	{
		locStepParticles.clear();
		dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_FinalParticles_Measured(Particle_t locStepInitialPID, Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();

	deque<const DParticleComboStep*> locParticleComboStepDeque;
	Get_ParticleComboSteps(locStepInitialPID, locParticleComboStepDeque);

	deque<const DKinematicData*> locStepParticles;
	for(size_t loc_i = 0; loc_i < locParticleComboStepDeque.size(); ++loc_i)
	{
		locStepParticles.clear();
		locParticleComboStepDeque[loc_i]->Get_FinalParticles_Measured(locPID, locStepParticles);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
}

void DParticleCombo::Get_DecayChainParticles_Measured(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const
{
	if((locStepIndex < 0) || (locStepIndex >= int(dParticleComboSteps.size())))
		return;
	const DParticleComboStep* locParticleComboStep = dParticleComboSteps[locStepIndex];

	if(locParticleComboStep->Get_InitialParticle_Measured() != NULL)
		locMeasuredParticles.push_back(locParticleComboStep->Get_InitialParticle_Measured());
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
			Get_DecayChainParticles_Measured(locParticleComboStep->Get_DecayStepIndex(loc_i), locMeasuredParticles);
		else if(locParticleComboStep->Is_FinalParticleDetected(loc_i))
			locMeasuredParticles.push_back(locParticleComboStep->Get_FinalParticle_Measured(loc_i));
	}
}


