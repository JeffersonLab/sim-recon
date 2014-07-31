#include "ANALYSIS/DReaction.h"

DReaction::DReaction(string locReactionName) : dReactionName(locReactionName)
{
	dKinFitType = d_NoFit;
	dChargedTrackFactoryTag = "";
	dNeutralShowerFactoryTag = "";

	dMinCombinedPIDFOM.first = false;
	dMinCombinedPIDFOM.second = 0.0;

	dMinCombinedTrackingFOM.first = false;
	dMinCombinedTrackingFOM.second = 0.0;

	dMinChargedPIDFOM.first = false;
	dMinChargedPIDFOM.second = 0.0;

	dMinPhotonPIDFOM.first = false;
	dMinPhotonPIDFOM.second = 5.73303E-7;

	dMinTrackingFOM.first = false;
	dMinTrackingFOM.second = 0.0027;

	dMaxPhotonRFDeltaT.first = false;
	dMaxPhotonRFDeltaT.second = 10.0*2.004; //10 RF buckets

	dMinProtonMomentum.first = true;
	dMinProtonMomentum.second = 0.25;

	dHasDetectorMatchFlag.first = false;
	dHasDetectorMatchFlag.second = false;

	dTTreeOutputFileName = "";
	dEnableTTreeOutputFlag = false;
	dMinThrownMatchFOMForROOT = -1.0; //always
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

string DReaction::Get_DecayChainFinalParticlesROOTNames(size_t locStepIndex, int locUpToStepIndex, deque<Particle_t> locUpThroughPIDs, bool locKinFitResultsFlag, bool locExpandDecayingParticlesFlag) const
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances)
	string locName = "";
	deque<Particle_t> locPIDs;
	dReactionSteps[locStepIndex]->Get_FinalParticleIDs(locPIDs);
	int locMissingParticleIndex = dReactionSteps[locStepIndex]->Get_MissingParticleIndex();
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		if(int(loc_j) == locMissingParticleIndex)
			continue; //exclude missing!

		if(int(locStepIndex) == locUpToStepIndex)
		{
			bool locPIDFoundFlag = false;
			for(deque<Particle_t>::iterator locIterator = locUpThroughPIDs.begin(); locIterator != locUpThroughPIDs.end(); ++locIterator)
			{
				if((*locIterator) != locPIDs[loc_j])
					continue;
				locUpThroughPIDs.erase(locIterator);
				locPIDFoundFlag = true;
				break;
			}
			if(!locPIDFoundFlag)
				continue; //skip it: don't want to include it
		}

		//find all future steps in which this pid is a parent
		int locDecayingStepIndex = -1;
		//expand if: not-kinfitting, or not a fixed mass
		if((!locKinFitResultsFlag) || (!IsFixedMass(locPIDs[loc_j])) || locExpandDecayingParticlesFlag)
		{
			for(size_t loc_k = locStepIndex + 1; loc_k < dReactionSteps.size(); ++loc_k)
			{
				if(dReactionSteps[loc_k]->Get_InitialParticleID() != locPIDs[loc_j])
					continue;
				locDecayingStepIndex = loc_k;
				break;
			}
		}

		if(locDecayingStepIndex == -1)
			locName += ParticleName_ROOT(locPIDs[loc_j]);
		else
			locName += Get_DecayChainFinalParticlesROOTNames(locDecayingStepIndex, locUpToStepIndex, locUpThroughPIDs, locKinFitResultsFlag, locExpandDecayingParticlesFlag);
	}
	return locName;
}

void DReaction::Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs, bool locIncludeDuplicatesFlag) const //independent of step
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

			if(!locIncludeDuplicatesFlag)
			{
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
			}			

			locDetectedPIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalChargedPIDs(deque<Particle_t>& locDetectedChargedPIDs, bool locIncludeDuplicatesFlag) const //independent of step
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

			if(!locIncludeDuplicatesFlag)
			{
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
			}

			locDetectedChargedPIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_FinalStatePIDs(deque<Particle_t>& locFinalStatePIDs, bool locIncludeDuplicatesFlag) const
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
			if(!locIncludeDuplicatesFlag)
			{
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
			}

			locFinalStatePIDs.push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs, bool locIncludeDuplicatesFlag) const
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

			if(!locIncludeDuplicatesFlag)
			{
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
			}

			locDetectedPIDs[loc_i].push_back(locPID);
		}
	}
}

void DReaction::Get_DetectedFinalChargedPIDs(deque<deque<Particle_t> >& locDetectedChargedPIDs, bool locIncludeDuplicatesFlag) const
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

			if(!locIncludeDuplicatesFlag)
			{
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
			}

			locDetectedChargedPIDs[loc_i].push_back(locPID);
		}
	}
}

