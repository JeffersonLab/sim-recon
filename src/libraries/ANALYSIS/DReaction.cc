#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisAction.h"

namespace DAnalysis
{

/************************************************************** DREACTION FUNCTIONS ***************************************************************/

vector<Particle_t> DReaction::Get_FinalPIDs(bool locIncludeMissingFlag, bool locIncludeDecayingFlag, Charge_t locCharge, bool locIncludeDuplicatesFlag) const
{
	vector<Particle_t> locFinalPIDs;
	for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
	{
		const DReactionStep* locStep = dReactionSteps[loc_i];
		for(size_t loc_j = 0; loc_j < locStep->Get_NumFinalPIDs(); ++loc_j)
		{
			if(!locIncludeDecayingFlag && (Get_DecayStepIndex(this, loc_i, loc_j) >= 0))
				continue;
			if(!locIncludeMissingFlag && (loc_j == locStep->Get_MissingParticleIndex()))
				continue;
			Particle_t locPID = locStep->Get_FinalPID(loc_j);
			if(Is_CorrectCharge(locPID, locCharge))
				locFinalPIDs.push_back(locPID);
		}
	}

	if(!locIncludeDuplicatesFlag)
	{
		std::sort(locFinalPIDs.begin(), locFinalPIDs.end()); //must sort first or else std::unique won't do what we want!
		locFinalPIDs.erase(std::unique(locFinalPIDs.begin(), locFinalPIDs.end()), locFinalPIDs.end());
	}
	return locFinalPIDs;
}

/************************************************************** NAMESPACE-SCOPE FUNCTIONS ***************************************************************/

pair<int, int> Get_InitialParticleDecayFromIndices(const DReaction* locReaction, int locStepIndex)
{
	//check to see how many initial-state particles with this PID type there are before now
	Particle_t locDecayingPID = locReaction->Get_ReactionStep(locStepIndex)->Get_InitialParticleID();
	size_t locPreviousPIDCount = 0;
	for(int loc_i = 0; loc_i < locStepIndex; ++loc_i)
	{
		if(locReaction->Get_ReactionStep(loc_i)->Get_InitialParticleID() == locDecayingPID)
			++locPreviousPIDCount;
	}

	//now, search through final-state PIDs until finding the (locPreviousPIDCount + 1)'th instance of this PID
	size_t locSearchPIDCount = 0;
	for(int loc_i = 0; loc_i < locStepIndex; ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			if(locReactionStep->Get_FinalParticleID(loc_j) != locDecayingPID)
				continue;
			++locSearchPIDCount;
			if(locSearchPIDCount <= locPreviousPIDCount)
				continue;
			return pair<int, int>(loc_i, loc_j);
		}
	}
	return pair<int, int>(-1, -1);
}

int Get_DecayStepIndex(const DReaction* locReaction, size_t locStepIndex, size_t locParticleIndex)
{
	//check if the input particle decays later in the reaction
	Particle_t locDecayingPID = locReaction->Get_ReactionStep(locStepIndex)->Get_FinalPID(locParticleIndex);

	//check to see how many final state particles with this pid type there are before now
	size_t locPreviousPIDCount = 0;
	for(size_t loc_i = 0; loc_i <= locStepIndex; ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			if((loc_i == locStepIndex) && (loc_j == locParticleIndex))
				break; //at the current particle: of the search
			if(locReactionStep->Get_FinalParticleID(loc_j) == locDecayingPID)
				++locPreviousPIDCount;
		}
	}

	//now, find the (locPreviousPIDCount + 1)'th time where this pid is a decay parent
	size_t locStepPIDCount = 0;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locReaction->Get_ReactionStep(loc_i)->Get_InitialParticleID() != locDecayingPID)
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

//CHANGE ME???
int Get_DefinedParticleStepIndex(const DReaction* locReaction) const
{
	//-1 if none //defined: missing or open-ended-decaying
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);

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

} //end namespace DAnalysis
