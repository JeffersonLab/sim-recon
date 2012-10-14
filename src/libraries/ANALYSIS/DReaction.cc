#include "ANALYSIS/DReaction.h"

const DReactionStep* DReaction::Get_ReactionStep(size_t locStepIndex) const
{
	if(locStepIndex >= dReactionSteps.size())
		return NULL;
	return dReactionSteps[locStepIndex];
}

string DReaction::Get_InitialParticlesROOTName(void) const
{
	if(dReactionSteps.empty())
		return (string());
	return dReactionSteps[0]->Get_InitialParticlesROOTName();
}

string DReaction::Get_DetectedParticlesROOTName(void) const
{
	string locDetectedParticlesROOTName;

	Particle_t locPID;
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			locDetectedParticlesROOTName += ParticleName_ROOT(locPID);
		}
	}
	return locDetectedParticlesROOTName;
}

void DReaction::Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<string>& locNames, bool locKinFitResultsFlag) const
{
	deque<deque<string> > locParticleNames;
	Get_DecayChainFinalParticlesROOTNames(locInitialPID, locParticleNames, locNames, locKinFitResultsFlag);
}

void DReaction::Get_DecayChainFinalParticlesROOTNames(Particle_t locInitialPID, deque<deque<string> >& locParticleNames, deque<string>& locNames, bool locKinFitResultsFlag) const
{
	locNames.clear();
	deque<Particle_t> locPIDs;
	locParticleNames.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		if(dReactionSteps[loc_i]->Get_InitialParticleID() != locInitialPID)
			continue;
		deque<deque<string> > locTempNames;
		Get_DecayChainFinalParticlesROOTNames(loc_i, locTempNames, locKinFitResultsFlag);
		locParticleNames.insert(locParticleNames.end(), locTempNames.begin(), locTempNames.end());
	}

	//eliminate duplicate names
		//particles could be a different order!!!
	deque<deque<string> >::iterator locIterator, locIterator2;
	deque<string>::iterator locIterator3, locIterator4;
	for(locIterator = locParticleNames.begin(); locIterator != locParticleNames.end(); ++locIterator)
	{
		deque<string> locName1 = *locIterator; //one group of particle names
		for(locIterator2 = locIterator + 1; locIterator2 != locParticleNames.end(); ++locIterator2)
		{
			deque<string> locName2 = *locIterator2; //another group of particle names
			if(locName1.size() != locName2.size())
				break; //not same size, clearly can't be the same

			//loop over the lists of particles, see if they're identical
			for(locIterator3 = locName1.begin(); locIterator3 != locName1.end(); ++locIterator3)
			{
				for(locIterator4 = locName2.begin(); locIterator4 != locName2.end(); ++locIterator4)
				{
					if((*locIterator3) == (*locIterator4))
					{
						locName2.erase(locIterator4); //particle name is identical, remove it from the list of remaining names
						break;
					}
				}
			}
			if(locName2.empty()) //all names removed means all names matched: duplicate
			{
				locIterator2 = locParticleNames.erase(locIterator2);
				--locIterator2;
			}
		}
	}

	//finally build the strings
	for(size_t loc_i = 0; loc_i < locParticleNames.size(); ++loc_i)
	{
		string locName;
		for(size_t loc_j = 0; loc_j < locParticleNames[loc_i].size(); ++loc_j)
			locName += locParticleNames[loc_i][loc_j];
		locNames.push_back(locName);
	}
}

void DReaction::Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, deque<deque<string> >& locNames, bool locKinFitResultsFlag) const
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances and excluded particles)
	deque<Particle_t> locPIDs;
	dReactionSteps[locStepIndex]->Get_FinalParticleIDs(locPIDs);
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		//find all future steps in which this pid is a parent
		deque<size_t> locDecayingStepIndices;
		for(size_t loc_k = locStepIndex + 1; loc_k < dReactionSteps.size(); ++loc_k)
		{
			if(dReactionSteps[loc_k]->Get_InitialParticleID() != locPIDs[loc_j])
				continue;
			//expand if: not-kinfitting, not a fixed mass, OR excluded from kinfit
			if((!locKinFitResultsFlag) || (!IsFixedMass(locPIDs[loc_j])) || Check_IfDecayingParticleExcludedFromP4KinFit(loc_k))
				locDecayingStepIndices.push_back(loc_k);
		}

		if(locDecayingStepIndices.empty())
		{
			if(locNames.empty())
				locNames.push_back(deque<string>());
			for(size_t loc_k = 0; loc_k < locNames.size(); ++loc_k)
				locNames[loc_k].push_back(ParticleName_ROOT(locPIDs[loc_j]));
			continue;
		}

		deque<deque<string> > locBaseNames = locNames;
		locNames.clear();
		for(size_t loc_k = 0; loc_k < locDecayingStepIndices.size(); ++loc_k)
		{
			deque<deque<string> > locTempNames = locBaseNames;
			Get_DecayChainFinalParticlesROOTNames(locDecayingStepIndices[loc_k], locTempNames, locKinFitResultsFlag);
			locNames.insert(locNames.end(), locTempNames.begin(), locTempNames.end());
		}
	}
}

bool DReaction::Check_IfDecayingParticleExcludedFromP4KinFit(size_t locStepIndex) const
{
	for(size_t loc_i = 0; loc_i < dDecayingParticlesExcludedFromP4Kinfit.size(); ++loc_i)
	{
		if(dDecayingParticlesExcludedFromP4Kinfit[loc_i] == locStepIndex)
			return true;
	}
	return false;
}

bool DReaction::Check_IsDecayingParticle(Particle_t locPID, size_t locSearchStartIndex) const
{
	//see if this pid is a parent in a future step
	for(size_t loc_k = locSearchStartIndex; loc_k < dReactionSteps.size(); ++loc_k)
	{
		if(dReactionSteps[loc_k]->Get_InitialParticleID() == locPID)
			return true;
	}
	return false;
}

void DReaction::Get_ReactionSteps(Particle_t locInitialPID, deque<const DReactionStep*>& locReactionSteps) const
{
	locReactionSteps.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		if(dReactionSteps[loc_i]->Get_InitialParticleID() == locInitialPID)
			locReactionSteps.push_back(dReactionSteps[loc_i]);
	}
}

void DReaction::Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs) const //independent of step
{
	Particle_t locPID;
	locDetectedPIDs.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//see if this PID is already stored
			bool locAlreadyHavePIDFlag = false;
			for(size_t loc_k = 0; loc_k < locDetectedPIDs.size(); ++loc_k)
			{
				if(locDetectedPIDs[loc_k] == locPID)
				{
					locAlreadyHavePIDFlag = true;
					break;
				}
			}
			if(locAlreadyHavePIDFlag)
				continue;

			locDetectedPIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalChargedPIDs(deque<Particle_t>& locDetectedChargedPIDs) const //independent of step
{
	Particle_t locPID;
	locDetectedChargedPIDs.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(ParticleCharge(locPID) == 0)
				continue; //neutral
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//see if this PID is already stored
			bool locAlreadyHavePIDFlag = false;
			for(size_t loc_k = 0; loc_k < locDetectedChargedPIDs.size(); ++loc_k)
			{
				if(locDetectedChargedPIDs[loc_k] == locPID)
				{
					locAlreadyHavePIDFlag = true;
					break;
				}
			}
			if(locAlreadyHavePIDFlag)
				continue;

			locDetectedChargedPIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_FinalStatePIDs(deque<Particle_t>& locFinalStatePIDs) const
{
	Particle_t locPID;
	locFinalStatePIDs.resize(0);
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//see if this PID is already stored
			bool locAlreadyHavePIDFlag = false;
			for(size_t loc_k = 0; loc_k < locFinalStatePIDs.size(); ++loc_k)
			{
				if(locFinalStatePIDs[loc_k] == locPID)
				{
					locAlreadyHavePIDFlag = true;
					break;
				}
			}
			if(locAlreadyHavePIDFlag)
				continue;

			locFinalStatePIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs) const
{
	Particle_t locPID;
	locDetectedPIDs.clear();
	locDetectedPIDs.resize(dReactionSteps.size());
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//see if this PID is already stored
			bool locAlreadyHavePIDFlag = false;
			for(size_t loc_k = 0; loc_k < locDetectedPIDs[loc_i].size(); ++loc_k)
			{
				if(locDetectedPIDs[loc_i][loc_k] == locPID)
				{
					locAlreadyHavePIDFlag = true;
					break;
				}
			}
			if(locAlreadyHavePIDFlag)
				continue;

			locDetectedPIDs[loc_i].push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalChargedPIDs(deque<deque<Particle_t> >& locDetectedChargedPIDs) const
{
	Particle_t locPID;
	locDetectedChargedPIDs.clear();
	locDetectedChargedPIDs.resize(dReactionSteps.size());
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(ParticleCharge(locPID) == 0)
				continue; //neutral
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			//see if this pid is a parent in a future step
			if(Check_IsDecayingParticle(locPID, loc_i + 1))
				continue;

			//see if this PID is already stored
			bool locAlreadyHavePIDFlag = false;
			for(size_t loc_k = 0; loc_k < locDetectedChargedPIDs[loc_i].size(); ++loc_k)
			{
				if(locDetectedChargedPIDs[loc_i][loc_k] == locPID)
				{
					locAlreadyHavePIDFlag = true;
					break;
				}
			}
			if(locAlreadyHavePIDFlag)
				continue;

			locDetectedChargedPIDs[loc_i].push_back(locPID);
		}
	}
}

DAnalysisAction* DReaction::Get_AnalysisAction(size_t locActionIndex) const
{
	if(locActionIndex >= dAnalysisActions.size())
		return NULL;
	return dAnalysisActions[locActionIndex];
}


