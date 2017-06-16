#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisAction.h"

namespace DAnalysis
{

/************************************************************** DREACTION FUNCTIONS ***************************************************************/

DReaction::~DReaction(void)
{
	//DO NOT DELETE REACTION STEPS: MIGHT BE SHARED BETWEEN DIFFERENT DREACTIONS
	for(const auto& locAction : dAnalysisActions)
		delete locAction;
}

vector<Particle_t> DReaction::Get_FinalPIDs(int locStepIndex, bool locIncludeMissingFlag, bool locIncludeDecayingFlag, Charge_t locCharge, bool locIncludeDuplicatesFlag) const
{
	//define the PID loop
	vector<Particle_t> locFinalPIDs;
	auto locPIDLoop = [&](const DReactionStep* locStep, size_t locLoopStepIndex) -> void
	{
		for(size_t locPIDIndex = 0; locPIDIndex < locStep->Get_NumFinalPIDs(); ++locPIDIndex)
		{
			if(!locIncludeDecayingFlag && (Get_DecayStepIndex(this, locLoopStepIndex, locPIDIndex) >= 0))
				continue;
			if(!locIncludeMissingFlag && (locPIDIndex == locStep->Get_MissingParticleIndex()))
				continue;
			Particle_t locPID = locStep->Get_FinalPID(locPIDIndex);
			if(Is_CorrectCharge(locPID, locCharge))
				locFinalPIDs.push_back(locPID);
		}
	};

	//execute the loop
	if(locStepIndex != -1)
		locPIDLoop(Get_ReactionStep(locStepIndex), locStepIndex);
	else
	{
		for(size_t loc_i = 0; loc_i < dReactionSteps.size(); ++loc_i)
			locPIDLoop(Get_ReactionStep(locStepIndex), loc_i);
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
	Particle_t locDecayingPID = locReaction->Get_ReactionStep(locStepIndex)->Get_InitialPID();

	auto locSteps = locReaction->Get_ReactionSteps();
	auto locPIDCounter = [=](const DReactionStep* locStep) -> bool {return (locStep->Get_InitialPID() == locDecayingPID);};
	size_t locPreviousPIDCount = std::count_if(locSteps.begin(), locSteps.begin() + locStepIndex, locPIDCounter);

	//now, search through final-state PIDs until finding the (locPreviousPIDCount + 1)'th instance of this PID
	size_t locSearchPIDCount = 0;
	for(size_t loc_i = 0; loc_i < locStepIndex; ++loc_i)
	{
		auto locFinalPIDs = locSteps[loc_i]->Get_FinalPIDs(false);

		auto locPIDCounter = [&](Particle_t locPID) -> bool
			{return (locPID != locDecayingPID) ? false : (++locSearchPIDCount > locPreviousPIDCount);};

		auto locIterator = std::find_if(locFinalPIDs.begin(), locFinalPIDs.end(), locPIDCounter);
		if(locIterator != locFinalPIDs.end())
			return pair<int, int>(loc_i, std::distance(locFinalPIDs.begin(), locIterator));
	}
	return make_pair(-1, -1);
}

size_t Get_ParticleInstanceIndex(const DReactionStep* locStep, size_t locParticleIndex)
{
	auto locFinalPIDs = locStep->Get_FinalPIDs(false);
	Particle_t locPID = locFinalPIDs[locParticleIndex];
	return std::count(locFinalPIDs.begin(), locFinalPIDs.begin() + locParticleIndex, locPID) - 1; //-1: index starting from 0
}

int Get_DecayStepIndex(const DReaction* locReaction, size_t locStepIndex, size_t locParticleIndex)
{
	//check if the input particle decays later in the reaction
	auto locSteps = locReaction->Get_ReactionSteps();
	Particle_t locDecayingPID = locSteps[locStepIndex]->Get_FinalPID(locParticleIndex);

	//check to see how many final state particles with this pid type there are before now
	size_t locPreviousPIDCount = 0;
	for(size_t loc_i = 0; loc_i <= locStepIndex; ++loc_i)
	{
		const DReactionStep* locStep = locSteps[loc_i];
		auto locFinalPIDs = locStep->Get_FinalPIDs();
		auto locEndIterator = (loc_i++ == locStepIndex) ? locFinalPIDs.begin() + locParticleIndex : locFinalPIDs.end();
		locPreviousPIDCount += std::count(locFinalPIDs.begin(), locEndIterator, locDecayingPID);
	}

	//now, find the (locPreviousPIDCount + 1)'th time where this pid is a decay parent
	size_t locStepPIDCount = 0;
	auto locStepCounter = [&](const DReactionStep* locStep) -> bool
		{return (locStep->Get_InitialPID() != locDecayingPID) ? false : (++locStepPIDCount > locPreviousPIDCount);};

	auto locIterator = std::find_if(locSteps.begin(), locSteps.end(), locStepCounter);
	return (locIterator == locSteps.end()) ? -1 : std::distance(locSteps.begin(), locIterator);
}

vector<Particle_t> Get_ChainPIDs(const DReaction* locReaction, size_t locStepIndex, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locExpandDecayingFlag)
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances)
	vector<Particle_t> locChainPIDs;
	auto locPIDs = locReaction->Get_ReactionStep(locStepIndex)->Get_FinalPIDs();
	for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
	{
		Particle_t locPID = locPIDs[loc_j];
		if(!locUpThroughPIDs.empty() && (int(locStepIndex) == locUpToStepIndex))
		{
			auto locIterator = std::find(locUpThroughPIDs.begin(), locUpThroughPIDs.end(), locPID);
			if(locIterator != locUpThroughPIDs.end())
				locUpThroughPIDs.erase(locIterator); //found, include it
			else
				continue; //skip it: don't want to include it
		}

		//see if decaying, and if so where, if decay expansion is requested (or mass isn't fixed)
		int locDecayingStepIndex = (!IsFixedMass(locPID) || locExpandDecayingFlag) ? Get_DecayStepIndex(locReaction, locStepIndex, loc_j) : -1;

		if(locDecayingStepIndex == -1)
			locChainPIDs.push_back(locPID);
		else
		{
			auto locDecayPIDs = Get_ChainPIDs(locDecayingStepIndex, locUpToStepIndex, locUpThroughPIDs, locExpandDecayingFlag);
			locChainPIDs.insert(locChainPIDs.end(), locDecayPIDs.begin(), locDecayPIDs.end());
		}
	}
	return locChainPIDs;
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
