#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisAction.h"

DReaction::DReaction(string locReactionName) : dReactionName(locReactionName)
{
	dKinFitType = d_NoFit;

	dMinProtonMomentum = pair<bool, double>(true, 0.25);
	dMinChargedPIDFOM = pair<bool, double>(false, 5.73303E-7);
	dMinPhotonPIDFOM = pair<bool, double>(false, 5.73303E-7);
	dMaxPhotonRFDeltaT = pair<bool, double>(false, 0.5*2.004);
	dMaxExtraGoodTracks = pair<bool, size_t>(false, 4);
	dMaxNumBeamPhotonsInBunch = pair<bool, size_t>(false, 0);

	dTTreeOutputFileName = "";
	dEnableTTreeOutputFlag = false;

	dEventStoreQuery = pair<string, string>("all", "");
}

DReaction::~DReaction(void)
{
	//DO NOT DELETE REACTION STEPS: MIGHT BE SHARED BETWEEN DIFFERENT DREACTIONS
	for(size_t loc_i = 0; loc_i < dAnalysisActions.size(); ++loc_i)
		delete dAnalysisActions[loc_i];
	for(size_t loc_i = 0; loc_i < dComboPreSelectionActions.size(); ++loc_i)
		delete dComboPreSelectionActions[loc_i];
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

void DReaction::Get_DetectedFinalPIDs(deque<Particle_t>& locDetectedPIDs, int locChargeFlag, bool locIncludeDuplicatesFlag) const //independent of step
{
	//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles

	locDetectedPIDs.clear();
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			Particle_t locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			int locCharge = ParticleCharge(locPID);
			if(((locChargeFlag == 1) && (locCharge == 0)) || ((locChargeFlag == 2) && (locCharge != 0)))
				continue; //wrong charge
			if(((locChargeFlag == 3) && (locCharge <= 0)) || ((locChargeFlag == 4) && (locCharge >= 0)))
				continue; //wrong charge

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

void DReaction::Get_DetectedFinalPIDs(deque<deque<Particle_t> >& locDetectedPIDs, int locChargeFlag, bool locIncludeDuplicatesFlag) const
{
	//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles

	locDetectedPIDs.clear();
	locDetectedPIDs.resize(dReactionSteps.size());
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < dReactionSteps[loc_i]->Get_NumFinalParticleIDs(); ++loc_j)
		{
			Particle_t locPID = dReactionSteps[loc_i]->Get_FinalParticleID(loc_j);
			if(dReactionSteps[loc_i]->Get_MissingParticleIndex() == int(loc_j))
				continue; //missing particle

			int locCharge = ParticleCharge(locPID);
			if(((locChargeFlag == 1) && (locCharge == 0)) || ((locChargeFlag == 2) && (locCharge != 0)))
				continue; //wrong charge
			if(((locChargeFlag == 3) && (locCharge <= 0)) || ((locChargeFlag == 4) && (locCharge >= 0)))
				continue; //wrong charge

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

int DReaction::Get_DecayStepIndex(int locStepIndex, int locParticleIndex) const
{
	//check if the ionput particle decays later in the reaction
	Particle_t locDecayingPID = Get_ReactionStep(locStepIndex)->Get_FinalParticleID(locParticleIndex);

	if((locDecayingPID == Gamma) || (locDecayingPID == Electron) || (locDecayingPID == Positron) || (locDecayingPID == Proton) || (locDecayingPID == AntiProton))
		return -1; //these particles don't decay: don't search!

	//check to see how many final state particles with this pid type there are before now
	size_t locPreviousPIDCount = 0;
	for(int loc_i = 0; loc_i <= locStepIndex; ++loc_i)
	{
		const DReactionStep* locReactionStep = Get_ReactionStep(loc_i);
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			if((loc_i == locStepIndex) && (int(loc_j) == locParticleIndex))
				break; //at the current particle: of the search
			if(locReactionStep->Get_FinalParticleID(loc_j) == locDecayingPID)
				++locPreviousPIDCount;
		}
	}

	//now, find the (locPreviousPIDCount + 1)'th time where this pid is a decay parent
	size_t locStepPIDCount = 0;
	for(size_t loc_i = 0; loc_i < Get_NumReactionSteps(); ++loc_i)
	{
		if(Get_ReactionStep(loc_i)->Get_InitialParticleID() != locDecayingPID)
			continue;
		++locStepPIDCount;
		if(locStepPIDCount <= locPreviousPIDCount)
			continue;
		//decays later in the reaction, at step index loc_i
		return loc_i;
	}

	// does not decay later in the reaction
	return -1;
}

int DReaction::Get_DefinedParticleStepIndex(void) const
{
	//-1 if none //defined: missing or open-ended-decaying
	for(size_t loc_i = 0; loc_i < Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = Get_ReactionStep(loc_i);

		//check for open-ended-decaying particle
		Particle_t locTargetPID = locReactionStep->Get_TargetParticleID();
		if((loc_i == 0) && (locTargetPID == Unknown))
			return loc_i;

		//check for missing particle
		Particle_t locMissingPID = Unknown;
		if(locReactionStep->Get_MissingPID(locMissingPID))
			return loc_i;
	}

	return -1;
}


