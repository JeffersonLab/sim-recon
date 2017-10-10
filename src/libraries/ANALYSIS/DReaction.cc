#include "ANALYSIS/DReaction.h"
#include "ANALYSIS/DAnalysisAction.h"
#include "ANALYSIS/DCutActions.h"

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
	vector<Particle_t> locFinalPIDs;
	for(size_t locLoopStepIndex = 0; locLoopStepIndex < dReactionSteps.size(); ++locLoopStepIndex)
	{
		if((locStepIndex != -1) && (int(locLoopStepIndex) != locStepIndex))
			continue;

		auto locStep = dReactionSteps[locLoopStepIndex];
		for(size_t locPIDIndex = 0; locPIDIndex < locStep->Get_NumFinalPIDs(); ++locPIDIndex)
		{
			if(!locIncludeDecayingFlag && (DAnalysis::Get_DecayStepIndex(this, locLoopStepIndex, locPIDIndex) >= 0))
				continue;
			if(!locIncludeMissingFlag && (int(locPIDIndex) == locStep->Get_MissingParticleIndex()))
				continue;
			Particle_t locPID = locStep->Get_FinalPID(locPIDIndex);
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

vector<Particle_t> DReaction::Get_MissingPIDs(int locStepIndex, Charge_t locCharge, bool locIncludeDuplicatesFlag) const
{
	vector<Particle_t> locFinalPIDs;
	for(size_t locLoopStepIndex = 0; locLoopStepIndex < dReactionSteps.size(); ++locLoopStepIndex)
	{
		if((locStepIndex != -1) && (int(locLoopStepIndex) != locStepIndex))
			continue;

		auto locStep = dReactionSteps[locLoopStepIndex];
		for(size_t locPIDIndex = 0; locPIDIndex < locStep->Get_NumFinalPIDs(); ++locPIDIndex)
		{
			if(int(locPIDIndex) != locStep->Get_MissingParticleIndex())
				continue;
			Particle_t locPID = locStep->Get_FinalPID(locPIDIndex);
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

void DReaction::Set_MaxPhotonRFDeltaT(double locMaxPhotonRFDeltaT)
{
	cout << "WARNING: USING DReaction::Set_MaxPhotonRFDeltaT() IS DEPRECATED. PLEASE SWITCH TO Set_NumPlusMinusRFBunches()." << endl;
	dMaxPhotonRFDeltaT = pair<bool, double>(true, locMaxPhotonRFDeltaT);
}

void DReaction::Add_ComboPreSelectionAction(DAnalysisAction* locAction)
{
	cout << "WARNING: USING DReaction::Add_ComboPreSelectionAction() IS DEPRECATED. PLEASE SWITCH TO Add_AnalysisAction()." << endl;
	Add_AnalysisAction(locAction);
}

void DReaction::Set_InvariantMassCut(Particle_t locStepInitialPID, double locMinInvariantMass, double locMaxInvariantMass)
{
	cout << "WARNING: USING DReaction::Set_InvariantMassCut() IS DEPRECATED. PLEASE SWITCH TO Add_AnalysisAction()." << endl;
	Add_AnalysisAction(new DCutAction_InvariantMass(this, locStepInitialPID, false, locMinInvariantMass, locMaxInvariantMass));
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
	for(int loc_i = 0; loc_i < locStepIndex; ++loc_i)
	{
		auto locFinalPIDs = locSteps[loc_i]->Get_FinalPIDs(false); //exclude missing

		auto locPIDCounter = [&](Particle_t locPID) -> bool
			{return (locPID != locDecayingPID) ? false : (++locSearchPIDCount > locPreviousPIDCount);};

		auto locIterator = std::find_if(locFinalPIDs.begin(), locFinalPIDs.end(), locPIDCounter);
		if(locIterator != locFinalPIDs.end())
			return pair<int, int>(loc_i, std::distance(locFinalPIDs.begin(), locIterator));
	}
	return std::make_pair(-1, -1);
}

size_t Get_ParticleInstanceIndex(const DReactionStep* locStep, size_t locParticleIndex)
{
	auto locFinalPIDs = locStep->Get_FinalPIDs();
	auto locPID = locFinalPIDs[locParticleIndex];
	size_t locInstanceIndex = 0;

	//see how many of same PID came before the particle index
	for(size_t loc_i = 0; loc_i < locParticleIndex; ++loc_i)
	{
		if(locStep->Get_MissingParticleIndex() == int(loc_i))
			continue;
		if(locFinalPIDs[loc_i] == locPID)
			++locInstanceIndex;
	}
	return locInstanceIndex;
}

int Get_DecayStepIndex(const DReaction* locReaction, size_t locStepIndex, size_t locParticleIndex)
{
	//check if the input particle decays later in the reaction
	auto locSteps = locReaction->Get_ReactionSteps();
	auto locDecayingPID = locSteps[locStepIndex]->Get_FinalPID(locParticleIndex);
	if(locSteps[locStepIndex]->Get_MissingParticleIndex() == int(locParticleIndex))
		return -1; //missing, does not decay

	//check to see how many final state particles with this pid type there are before now
	size_t locPreviousPIDCount = 0;
	for(size_t loc_i = 0; loc_i <= locStepIndex; ++loc_i)
	{
		auto locStep = locSteps[loc_i];
		auto locFinalPIDs = locStep->Get_FinalPIDs(false); //exclude missing
		auto locEndIterator = (loc_i == locStepIndex) ? locFinalPIDs.begin() + locParticleIndex : locFinalPIDs.end();
		locPreviousPIDCount += std::count(locFinalPIDs.begin(), locEndIterator, locDecayingPID);
	}

	//now, find the (locPreviousPIDCount + 1)'th time where this pid is a decay parent
	size_t locStepPIDCount = 0;
	auto StepCounter = [&](const DReactionStep* locStep) -> bool
		{return (locStep->Get_InitialPID() != locDecayingPID) ? false : (++locStepPIDCount > locPreviousPIDCount);};

	auto locIterator = std::find_if(std::next(locSteps.begin()), locSteps.end(), StepCounter);
	return (locIterator == locSteps.end()) ? -1 : std::distance(locSteps.begin(), locIterator);
}

vector<Particle_t> Get_ChainPIDs(const DReaction* locReaction, size_t locStepIndex, int locUpToStepIndex, vector<Particle_t> locUpThroughPIDs, bool locExpandDecayingFlag, bool locExcludeMissingFlag)
{
	//if locKinFitResultsFlag = true: don't expand decaying particles (through decay chain) that were included in the kinfit (still expand resonances)
	vector<Particle_t> locChainPIDs;
	auto locPIDs = locReaction->Get_ReactionStep(locStepIndex)->Get_FinalPIDs(!locExcludeMissingFlag);
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
			auto locDecayPIDs = Get_ChainPIDs(locReaction, locDecayingStepIndex, locUpToStepIndex, locUpThroughPIDs, locExpandDecayingFlag, locExcludeMissingFlag);
			locChainPIDs.insert(locChainPIDs.end(), locDecayPIDs.begin(), locDecayPIDs.end());
		}
	}
	return locChainPIDs;
}

vector<size_t> Get_DefinedParticleStepIndex(const DReaction* locReaction)
{
	//-1 if none, -2 if more than 1 //defined: missing or open-ended-decaying
	vector<size_t> locDefinedParticleStepIndices;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		auto locReactionStep = locReaction->Get_ReactionStep(loc_i);

		//check for open-ended-decaying particle
		if((loc_i == 0) && (locReactionStep->Get_TargetPID() == Unknown))
		{
			locDefinedParticleStepIndices.push_back(loc_i);
			continue;
		}

		//check for missing particle
		if(locReactionStep->Get_IsInclusiveFlag() || (locReactionStep->Get_MissingPID() != Unknown))
		{
			locDefinedParticleStepIndices.push_back(loc_i);
			continue;
		}
	}

	return locDefinedParticleStepIndices;
}

vector<const DReaction*> Get_Reactions(JEventLoop* locEventLoop)
{
	// Get DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	vector<const DReaction*> locReactions;
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>*>(locFactories[loc_i]);
		if(locFactory == nullptr)
			continue;
		if(string(locFactory->Tag()) == "Thrown")
			continue;

		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		locReactions.insert(locReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}
	return locReactions;
}

} //end namespace DAnalysis
