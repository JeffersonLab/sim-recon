#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DAnalysisUtilities.h"
#include "ANALYSIS/DParticleComboCreator.h"

DAnalysisUtilities::DAnalysisUtilities(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dPIDAlgorithm);

	dTargetZCenter = 65.0;
	// Get Target parameters from XML
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber()) : NULL;
	locGeometry->GetTargetZ(dTargetZCenter);

	//Get magnetic field map
	dMagneticFieldMap = locApplication->GetBfield(locEventLoop->GetJEvent().GetRunNumber());

	//For "Unused" tracks/showers
	//BEWARE: IF THIS IS CHANGED, CHANGE IN THE BLUEPRINT FACTORY AND THE EVENT WRITER ALSO!!
	dTrackSelectionTag = "PreSelect";
	dShowerSelectionTag = "PreSelect";
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);
}

bool DAnalysisUtilities::Check_IsBDTSignalEvent(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag, bool locIncludeDecayingToReactionFlag) const
{
#ifdef VTRACE
	VT_TRACER("DAnalysisUtilities::Check_IsBDTSignalEvent()");
#endif

	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//if locIncludeDecayingToReactionFlag = true, will test whether the thrown reaction could decay to the DReaction
		//Note that resonances, phi's, and omega's are automatically decayed
			//e.g. if DReaction or thrown is g, p -> pi+, pi-, omega, p; will instead treat it as g, p -> 2pi+, 2pi-, pi0, p (or whatever the omega decay products are)
	//e.g. g, p -> pi+, pi0, K0, Lambda can decay to g, p -> 2pi+, 2pi-, pi0, p
		//if locIncludeDecayingToReactionFlag = true, then it would be included as "Signal," if false, then background
		//locIncludeDecayingToReactionFlag should be true UNLESS you are explicitly checking all possible reactions that could decay to your channel in your BDT
			//e.g. could kinfit to g, p -> pi+, pi0, K0, Lambda and include it as a BDT variable

	if(dParticleComboCreator == nullptr) //Can't create in constructor: infinite recursion
		dParticleComboCreator = new DParticleComboCreator(locEventLoop, nullptr, nullptr, nullptr);
	DReaction_factory_Thrown* dThrownReactionFactory = static_cast<DReaction_factory_Thrown*>(locEventLoop->GetFactory("DReaction", "Thrown"));

	vector<const DReaction*> locThrownReactions;
	locEventLoop->Get(locThrownReactions, "Thrown");
	if(locThrownReactions.empty())
		return false;
	auto locActualThrownReaction = locThrownReactions[0];

	//Replace omega & phi in DReaction with their decay products
	size_t locStepIndex = 0;
	vector<DReactionStep> locNewReactionSteps;
	const DReaction* locCurrentReaction = locReaction;
	DReaction locNewReaction("Interim"); //if needed

	do
	{
		const DReactionStep* locReactionStep = locCurrentReaction->Get_ReactionStep(locStepIndex);
		Particle_t locInitialPID = locReactionStep->Get_InitialPID();
		if((locInitialPID != omega) && (locInitialPID != phiMeson))
		{
			++locStepIndex;
			continue;
		}

		//is decaying phi or omega, replace it with its decay products
		auto locDecayProducts = locReactionStep->Get_FinalPIDs();

		//find the production step
		for(size_t loc_i = 0; loc_i < locStepIndex; ++loc_i)
		{
			const DReactionStep* locProductionStep = locCurrentReaction->Get_ReactionStep(loc_i);
			auto locPIDs = locProductionStep->Get_FinalPIDs();

			//search for the decaying particle. when found, replace it with its decay products
			bool locFoundFlag = false;
			for(size_t loc_j = 0; loc_j < locPIDs.size(); ++loc_j)
			{
				Particle_t locFinalStatePID = locPIDs[loc_j];
				if(locFinalStatePID != locInitialPID)
					continue;
				locPIDs.erase(locPIDs.begin() + loc_j);
				locPIDs.insert(locPIDs.begin(), locDecayProducts.begin(), locDecayProducts.end());
				locFoundFlag = true;
				break;
			}
			if(!locFoundFlag)
				continue;

			//make a new reaction step
			DReactionStep locNewReactionStep;
			locNewReactionStep.Set_InitialParticleID(locProductionStep->Get_InitialPID());
			locNewReactionStep.Set_TargetParticleID(locProductionStep->Get_TargetPID());
			int locMissingParticleIndex = locProductionStep->Get_MissingParticleIndex();
			for(int loc_j = 0; loc_j < int(locPIDs.size()); ++loc_j)
				locNewReactionStep.Add_FinalParticleID(locPIDs[loc_j], (loc_j == locMissingParticleIndex));
			locNewReactionSteps.push_back(locNewReactionStep);

			//make a new reaction
			DReaction locReactionBuild("Build");
			locReactionBuild.Clear_ReactionSteps();
			for(size_t loc_j = 0; loc_j < locCurrentReaction->Get_NumReactionSteps(); ++loc_j)
			{
				if(loc_j == loc_i)
					locReactionBuild.Add_ReactionStep(&locNewReactionSteps.back());
				else if(loc_j == locStepIndex)
					continue;
				else
					locReactionBuild.Add_ReactionStep(locCurrentReaction->Get_ReactionStep(loc_j));
			}
			locNewReaction = locReactionBuild;
			locCurrentReaction = &locNewReaction;
			break;
		}
	}
	while(locStepIndex < locCurrentReaction->Get_NumReactionSteps());

	//Get thrown steps
	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	Get_ThrownParticleSteps(locEventLoop, locThrownSteps);

	if(locThrownSteps.size() == 1) //nothing to replace: it either works or it doesn't
		return Check_ThrownsMatchReaction(locEventLoop, locCurrentReaction, locExclusiveMatchFlag);

	//build maps of counts of the thrown & DReaction decaying particles
	map<Particle_t, size_t> locNumDecayingParticles_Thrown;
	for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
	{
		if(locThrownSteps[loc_i].first == NULL)
			continue; //beam particle
		Particle_t locPID = locThrownSteps[loc_i].first->PID();
		if(locNumDecayingParticles_Thrown.find(locPID) == locNumDecayingParticles_Thrown.end())
			locNumDecayingParticles_Thrown[locPID] = 1;
		else
			++locNumDecayingParticles_Thrown[locPID];
	}

	map<Particle_t, size_t> locNumDecayingParticles_Reaction;
	for(size_t loc_i = 0; loc_i < locCurrentReaction->Get_NumReactionSteps(); ++loc_i)
	{
		Particle_t locPID = locCurrentReaction->Get_ReactionStep(loc_i)->Get_InitialPID();
		if(locPID == Gamma)
			continue; //beam particle
		if(locNumDecayingParticles_Reaction.find(locPID) == locNumDecayingParticles_Reaction.end())
			locNumDecayingParticles_Reaction[locPID] = 1;
		else
			++locNumDecayingParticles_Reaction[locPID];
	}

	//if there is a particle type in the DReaction that is not present in the thrown, it's gonna fail no matter what: check now
	map<Particle_t, size_t>::iterator locIterator = locNumDecayingParticles_Reaction.begin();
	for(; locIterator != locNumDecayingParticles_Reaction.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		if(locNumDecayingParticles_Thrown.find(locPID) == locNumDecayingParticles_Thrown.end())
			return false;
	}

	//loop through thrown steps: replace phi's, omega's with their decay products
		//also replace particles not in the DReaction with their decay products IF locIncludeDecayingToReactionFlag = true
	locStepIndex = 1;
	do
	{
		pair<const DMCThrown*, deque<const DMCThrown*> > locStepPair = locThrownSteps[locStepIndex];

		//check to see if the thrown decaying particle is present in the dreaction
		const DMCThrown* locThrownParent = locStepPair.first;
		Particle_t locInitialPID = locThrownParent->PID();
		bool locInitialPIDFoundFlag = (locNumDecayingParticles_Reaction.find(locInitialPID) != locNumDecayingParticles_Reaction.end());

		//if it was not found, the only way it can be a match is if the decay products ARE in the reaction step: decay it in place
			//also: omega and phi have non-negligible width: cannot constrain their peaks in kinfit or BDT: replace with their decay products
		if(((!locInitialPIDFoundFlag) && locIncludeDecayingToReactionFlag) || (locInitialPID == omega) || (locInitialPID == phiMeson))
		{
			Replace_DecayingParticleWithProducts(locThrownSteps, locStepIndex);
			if(locNumDecayingParticles_Thrown[locInitialPID] == 1)
				locNumDecayingParticles_Thrown.erase(locInitialPID);
			else
				--locNumDecayingParticles_Thrown[locInitialPID];
		}
		else
			++locStepIndex;
	}
	while(locStepIndex < locThrownSteps.size());

	if(!locIncludeDecayingToReactionFlag)
	{
		//don't try decaying thrown particles: compare as-is
		DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locThrownSteps);
		auto locThrownCombo = dParticleComboCreator->Build_ThrownCombo(locEventLoop, locThrownReaction, locThrownSteps);
		bool locCheckResult = Check_ThrownsMatchReaction(locActualThrownReaction, locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

		dParticleComboCreator->Reset();
		dThrownReactionFactory->Recycle_Reaction(locThrownReaction);
		return locCheckResult;
	}

	//at this point, #-steps are final. If exclusive, bail if they don't match
	if(locExclusiveMatchFlag)
	{
		if(locReaction->Get_NumReactionSteps() != locThrownSteps.size())
			return false;
	}

	//ok, if there are still an unequal # of parents for a given PID between thrown & DReaction (i.e. both #'s are non-zero):
		//it's too confusing to figure out which should be replaced by their decay products and which shouldn't
		//so, try all possibilities, and see if any of them match. 

	//build PIDs-to-replace vectors //easiest to manipulate a 1D vector
	vector<Particle_t> locPIDVector;
	vector<int> locResumeAtIndex;
	for(locIterator = locNumDecayingParticles_Reaction.begin(); locIterator != locNumDecayingParticles_Reaction.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		size_t locNumThrown = locNumDecayingParticles_Thrown[locPID];
		if(locNumThrown < locIterator->second)
			return false; //thrown doesn't match
		if(locNumThrown == locIterator->second)
			continue; //all is well
		for(size_t loc_i = 0; loc_i < locNumThrown - locIterator->second; ++loc_i)
		{
			locPIDVector.push_back(locPID);
			locResumeAtIndex.push_back(0);
		}
	}

	//if no additional replacements to make: check it
	if(locPIDVector.empty())
	{
		DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locThrownSteps);
		auto locThrownCombo = dParticleComboCreator->Build_ThrownCombo(locEventLoop, locThrownReaction, locThrownSteps);
		bool locCheckResult = Check_ThrownsMatchReaction(locActualThrownReaction, locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

		dParticleComboCreator->Reset();
		dThrownReactionFactory->Recycle_Reaction(locThrownReaction);
		return locCheckResult;
	}

	//loop through, making all combos of particles, (almost) just like in DParticleComboBlueprint_factory
		//unlike that factory, the below may try a given combo twice: too hard to code against it though 
			//for a given PID, hard to compare current locResumeAtIndex to previous ones, since the thrown steps are continually replaced
	deque<int> locDecayReplacementIndices;
	int locParticleIndex = 0;
	//below: contains "locThrownSteps" after each replacement
		//1st deque index is replacement index, 
	deque<deque<pair<const DMCThrown*, deque<const DMCThrown*> > > > locReplacementThrownSteps(1, locThrownSteps);
	do
	{
		if(locParticleIndex == -1)
			break; //no success

		deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locCurrentThrownSteps = locReplacementThrownSteps.back();
		if(locParticleIndex == int(locPIDVector.size()))
		{
			//combo defined: try it
			DReaction* locThrownReaction = dThrownReactionFactory->Build_ThrownReaction(locEventLoop, locCurrentThrownSteps);
			auto locThrownCombo = dParticleComboCreator->Build_ThrownCombo(locEventLoop, locThrownReaction, locCurrentThrownSteps);
			bool locCheckResult = Check_ThrownsMatchReaction(locActualThrownReaction, locThrownCombo, locCurrentReaction, locExclusiveMatchFlag);

			dParticleComboCreator->Reset();
			dThrownReactionFactory->Recycle_Reaction(locThrownReaction);

			if(locCheckResult)
				return true; //it worked!
			locReplacementThrownSteps.pop_back();
			--locParticleIndex;
			continue;
		}

		Particle_t locPID = locPIDVector[locParticleIndex];
		//find the next instance of this step & replace it
		bool locFoundFlag = false;
		for(size_t loc_i = locResumeAtIndex[locParticleIndex]; loc_i < locCurrentThrownSteps.size(); ++loc_i)
		{
			if(locCurrentThrownSteps[loc_i].first == NULL)
				continue;
			Particle_t locInitialPID = locCurrentThrownSteps[loc_i].first->PID();
			if(locInitialPID != locPID)
				continue;
			locFoundFlag = true;

			Replace_DecayingParticleWithProducts(locCurrentThrownSteps, loc_i);
			locReplacementThrownSteps.push_back(locCurrentThrownSteps);
			locResumeAtIndex[locParticleIndex] = loc_i + 1;

			break;
		}

		if(locFoundFlag)
			++locParticleIndex;
		else
		{
			locResumeAtIndex[locParticleIndex] = 0;
			locReplacementThrownSteps.pop_back();
			--locParticleIndex;
		}
	}
	while(true);

	return false;
}

void DAnalysisUtilities::Replace_DecayingParticleWithProducts(deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps, size_t locStepIndex) const
{
	if(locStepIndex >= locThrownSteps.size())
		return;

	//find the step where this particle is produced at
	int locInitialParticleDecayFromStepIndex = -1;
	const DMCThrown* locThrownParent = locThrownSteps[locStepIndex].first;
	if(locThrownParent == NULL)
		return;

	for(size_t loc_j = 0; loc_j < locStepIndex; ++loc_j)
	{
		for(size_t loc_k = 0; loc_k < locThrownSteps[loc_j].second.size(); ++loc_k)
		{
			if(locThrownParent != locThrownSteps[loc_j].second[loc_k])
				continue;
			locInitialParticleDecayFromStepIndex = loc_j;
			break;
		}
		if(locInitialParticleDecayFromStepIndex != -1)
			break;
	}

	if(locInitialParticleDecayFromStepIndex == -1)
	{
		cout << "ERROR: SOMETHING IS WRONG WITH DAnalysisUtilities::Replace_DecayingParticleWithProducts(). ABORTING." << endl;
		abort();
	}

	//insert the decay products where the production occured
	deque<const DMCThrown*>& locProductionStepFinalParticles = locThrownSteps[locInitialParticleDecayFromStepIndex].second;
	deque<const DMCThrown*>& locDecayStepFinalParticles = locThrownSteps[locStepIndex].second;
	for(size_t loc_j = 0; loc_j < locProductionStepFinalParticles.size(); ++loc_j)
	{
		if(locProductionStepFinalParticles[loc_j] != locThrownParent)
			continue;
		locProductionStepFinalParticles.erase(locProductionStepFinalParticles.begin() + loc_j);
		locProductionStepFinalParticles.insert(locProductionStepFinalParticles.end(), locDecayStepFinalParticles.begin(), locDecayStepFinalParticles.end());
		break;
	}

	locThrownSteps.erase(locThrownSteps.begin() + locStepIndex);
}

bool DAnalysisUtilities::Check_ThrownsMatchReaction(JEventLoop* locEventLoop, const DReaction* locReaction, bool locExclusiveMatchFlag) const
{
	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//note, if you decay a final state particle (e.g. k+, pi+) in your input DReaction*, a match will NOT be found: the thrown reaction/combo is truncated
	//if locExclusiveMatchFlag = false, then allow the input DReaction to be a subset of the thrown
	if(dParticleComboCreator == nullptr) //Can't create in constructor: infinite recursion
		dParticleComboCreator = new DParticleComboCreator(locEventLoop, nullptr, nullptr, nullptr);
	auto locThrownCombo = dParticleComboCreator->Build_ThrownCombo(locEventLoop);

	vector<const DReaction*> locReactions;
	locEventLoop->Get(locReactions, "Thrown");
	if(locReactions.empty())
		return false;

	auto locResult = Check_ThrownsMatchReaction(locReactions[0], locThrownCombo, locReaction, locExclusiveMatchFlag);
	dParticleComboCreator->Reset();
	return locResult;
}

bool DAnalysisUtilities::Check_ThrownsMatchReaction(const DReaction* locThrownReaction, const DParticleCombo* locThrownCombo, const DReaction* locReaction, bool locExclusiveMatchFlag) const
{
#ifdef VTRACE
	VT_TRACER("DAnalysisUtilities::Check_ThrownsMatchReaction()");
#endif

	//IF DREACTION HAS A MISSING UNKNOWN PARTICLE, MUST USE locExclusiveMatchFlag = false

	//note, if you decay a final state particle (e.g. k+, pi+) in your input DReaction*, a match will NOT be found: the thrown reaction/combo is truncated
	//if locExclusiveMatchFlag = false, then allow the input DReaction to be a subset of the thrown

	if(locThrownCombo == NULL)
		return false;
	if(locExclusiveMatchFlag)
	{
		if(locReaction->Get_NumReactionSteps() != locThrownCombo->Get_NumParticleComboSteps())
			return false;
	}

	//build map of InitialParticleDecayFromStepIndex's for input reaction: assume that if it matters, the user wanted them in the input order
	map<size_t, int> locReactionInitialParticleDecayFromStepIndexMap; //first is step index (of locReaction) where particle is initial, second is where is final state particle
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		if(loc_i == 0)
		{
			if(locReactionStep->Get_InitialPID() == Gamma)
				locReactionInitialParticleDecayFromStepIndexMap[0] = -1;
		}
		//loop over final state particles, and if decaying, find the step they are a parent in
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalPIDs(); ++loc_j)
		{
			Particle_t locFinalStatePID = locReactionStep->Get_FinalPID(loc_j);
			//see if this pid is a parent in a future step
			for(size_t loc_k = loc_i; loc_k < locReaction->Get_NumReactionSteps(); ++loc_k)
			{
				if(locReaction->Get_ReactionStep(loc_k)->Get_InitialPID() != locFinalStatePID)
					continue;
				if(locReactionInitialParticleDecayFromStepIndexMap.find(loc_k) != locReactionInitialParticleDecayFromStepIndexMap.end())
					continue; //this step already accounted for
				locReactionInitialParticleDecayFromStepIndexMap[loc_k] = loc_i;
				break;
			}
		}
	}

	//since throwns are more organized, loop through those and back-check to dreaction
	set<size_t> locMatchedInputStepIndices; //step indices in input locReaction that are already matched-to
	map<int, int> locStepMatching; //locReaction map to thrown step map
	map<int, int> locReverseStepMatching; //thrown map to locReaction step map
	for(size_t loc_i = 0; loc_i < locThrownCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		int locInitialParticleDecayFromStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(locThrownReaction, loc_i).first;

		//find where to start the search for this step in locReaction
		size_t locStartSearchIndex = 0; //locReaction could be a subset of the total; start at the beginning unless ...:
		if(locReverseStepMatching.find(locInitialParticleDecayFromStepIndex) != locReverseStepMatching.end())
			locStartSearchIndex = locReverseStepMatching[locInitialParticleDecayFromStepIndex] + 1; //parent step was matched, start search for this one after it

		//loop through locReaction and try to find this thrown step
		bool locMatchFoundFlag = false;
		int locPossibleMatchIndex = -1;
		for(size_t loc_j = locStartSearchIndex; loc_j < locReaction->Get_NumReactionSteps(); ++loc_j)
		{
			const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_j);
			if(locMatchedInputStepIndices.find(loc_j) != locMatchedInputStepIndices.end())
				continue; //this step was already accounted for

			//when not exact match: allow user step to have a missing unknown particle
			if(!DAnalysis::Check_ChannelEquality(locThrownReaction->Get_ReactionStep(loc_i), locReactionStep, false, !locExclusiveMatchFlag))
				continue; //particles aren't the same

			//ok, now check to make sure that the parent particle in this step was produced the same way in both thrown & locReaction
			if(locReactionInitialParticleDecayFromStepIndexMap.find(loc_j) == locReactionInitialParticleDecayFromStepIndexMap.end())
			{
				//this (loc_j) parent particle's production step wasn't listed in the locReaction: locReaction is probably a subset of the total
				//a match is possible but not certain: (e.g. locReaction is pi0, eta -> pi0, pi0 and this step (loc_i) is a pi0)
				locPossibleMatchIndex = loc_j;
				continue; //keep searching in case a different one should be used instead //will use this if another isn't found
			}

			int locReactionInitialParticleDecayFromStepIndex = locReactionInitialParticleDecayFromStepIndexMap[loc_j];
			if(locReactionInitialParticleDecayFromStepIndex != -1) //locReaction is not beam particle
			{
				if(locStepMatching.find(locReactionInitialParticleDecayFromStepIndex) == locStepMatching.end())
					continue; //the step in locReaction where this (loc_j) parent particle was produced was not mapped to the thrown steps yet: this is not the step we want
				int locReactionInitialParticleDecayFromStepIndexMappedBackToThrown = locStepMatching[locReactionInitialParticleDecayFromStepIndex];
				if(locInitialParticleDecayFromStepIndex != locReactionInitialParticleDecayFromStepIndexMappedBackToThrown)
					continue; //the decaying parent particle in this step (loc_j) comes from a different step in thrown (locInitialParticleDecayFromStepIndex)/reaction: continue
			}
			else if(locInitialParticleDecayFromStepIndex != -1)
				continue; //reaction is beam but thrown is not beam

			//finally, a match! register it
			locMatchedInputStepIndices.insert(loc_j);
			locStepMatching[loc_j] = loc_i;
			locReverseStepMatching[loc_i] = loc_j;
			locMatchFoundFlag = true;
			break;
		}

		if((!locMatchFoundFlag) && (locPossibleMatchIndex != -1))
		{
			//need to use the possible match
			locMatchedInputStepIndices.insert(locPossibleMatchIndex);
			locStepMatching[locPossibleMatchIndex] = loc_i;
			locReverseStepMatching[loc_i] = locPossibleMatchIndex;
			locMatchFoundFlag = true;
		}
		if(locExclusiveMatchFlag && (!locMatchFoundFlag))
			return false; //needed an exact match and it wasn't found: bail
	}
	if(locExclusiveMatchFlag)
		return true;

	//locReaction could be a subset of thrown: check if all locReaction steps found
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locMatchedInputStepIndices.find(loc_i) == locMatchedInputStepIndices.end())
			return false; //one of the input steps wasn't matched: abort!
	}

	return true;
}

void DAnalysisUtilities::Get_UnusedChargedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DChargedTrack*>& locUnusedChargedTracks) const
{
	locUnusedChargedTracks.clear();
	locEventLoop->Get(locUnusedChargedTracks, "Combo");
	std::sort(locUnusedChargedTracks.begin(), locUnusedChargedTracks.end());

	auto locChargedSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Charged);
	for(size_t loc_i = 0; loc_i < locChargedSourceObjects.size(); ++loc_i)
	{
		for(auto locIterator = locUnusedChargedTracks.begin(); locIterator != locUnusedChargedTracks.end();)
		{
			if(locChargedSourceObjects[loc_i] != *locIterator)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedChargedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedTimeBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackTimeBased*>& locUnusedTimeBasedTracks) const
{
	locUnusedTimeBasedTracks.clear();
	locEventLoop->Get(locUnusedTimeBasedTracks);

	vector<const DTrackTimeBased*> locComboTimeBasedTracks;
	locEventLoop->Get(locComboTimeBasedTracks, "Combo");
	locUnusedTimeBasedTracks.insert(locUnusedTimeBasedTracks.end(), locComboTimeBasedTracks.begin(), locComboTimeBasedTracks.end());

	auto locChargedSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Charged);
	for(size_t loc_i = 0; loc_i < locChargedSourceObjects.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		auto locChargedTrack = static_cast<const DChargedTrack*>(locChargedSourceObjects[loc_i]);
		for(auto locIterator = locUnusedTimeBasedTracks.begin(); locIterator != locUnusedTimeBasedTracks.end();)
		{
			if(locChargedTrack->candidateid != (*locIterator)->candidateid)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedTimeBasedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedWireBasedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackWireBased*>& locUnusedWireBasedTracks) const
{
	locUnusedWireBasedTracks.clear();
	locEventLoop->Get(locUnusedWireBasedTracks);

	auto locChargedSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Charged);
	for(size_t loc_i = 0; loc_i < locChargedSourceObjects.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		auto locChargedTrack = static_cast<const DChargedTrack*>(locChargedSourceObjects[loc_i]);
		for(auto locIterator = locUnusedWireBasedTracks.begin(); locIterator != locUnusedWireBasedTracks.end();)
		{
			if(locChargedTrack->candidateid != (*locIterator)->candidateid)
				++locIterator; //not-used (yet)
			else
				locIterator = locUnusedWireBasedTracks.erase(locIterator); //used
		}
	}
}

void DAnalysisUtilities::Get_UnusedTrackCandidates(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DTrackCandidate*>& locUnusedTrackCandidates) const
{
	locUnusedTrackCandidates.clear();
	locEventLoop->Get(locUnusedTrackCandidates);

	auto locChargedSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Charged);
	set<unsigned int> locUsedCandidateIndices;
	for(size_t loc_i = 0; loc_i < locChargedSourceObjects.size(); ++loc_i)
	{
		//only need the candidate id: same for all hypotheses for a given track
		auto locChargedTrack = static_cast<const DChargedTrack*>(locChargedSourceObjects[loc_i]);
		locUsedCandidateIndices.insert(locChargedTrack->candidateid - 1); //id = index + 1
	}

	for(int loc_i = locUnusedTrackCandidates.size() - 1; loc_i >= 0; --loc_i)
	{
		if(locUsedCandidateIndices.find(loc_i) != locUsedCandidateIndices.end())
			locUnusedTrackCandidates.erase(locUnusedTrackCandidates.begin() + loc_i);
	}
}

void DAnalysisUtilities::Get_UnusedNeutralShowers(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralShower*>& locUnusedNeutralShowers) const
{
	locUnusedNeutralShowers.clear();
	locEventLoop->Get(locUnusedNeutralShowers, dShowerSelectionTag.c_str());

	auto locNeutralSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Neutral);
	for(size_t loc_i = 0; loc_i < locNeutralSourceObjects.size(); ++loc_i)
	{
		auto locNeutralShower = static_cast<const DNeutralShower*>(locNeutralSourceObjects[loc_i]);
		for(auto locIterator = locUnusedNeutralShowers.begin(); locIterator != locUnusedNeutralShowers.end();)
		{
			if(locNeutralShower != *locIterator)
				++locIterator;
			else
				locIterator = locUnusedNeutralShowers.erase(locIterator);
		}
	}
}

void DAnalysisUtilities::Get_UnusedNeutralParticles(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, vector<const DNeutralParticle*>& locUnusedNeutralParticles) const
{
	locUnusedNeutralParticles.clear();
	locEventLoop->Get(locUnusedNeutralParticles, dShowerSelectionTag.c_str());

	auto locNeutralSourceObjects = locParticleCombo->Get_FinalParticle_SourceObjects(d_Neutral);
	for(size_t loc_i = 0; loc_i < locNeutralSourceObjects.size(); ++loc_i)
	{
		auto locNeutralShower = static_cast<const DNeutralShower*>(locNeutralSourceObjects[loc_i]);
		for(auto locIterator = locUnusedNeutralParticles.begin(); locIterator != locUnusedNeutralParticles.end();)
		{
			if(locNeutralShower != (*locIterator)->dNeutralShower)
				++locIterator;
			else
				locIterator = locUnusedNeutralParticles.erase(locIterator);
		}
	}
}

void DAnalysisUtilities::Get_ThrownParticleSteps(JEventLoop* locEventLoop, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps) const
{
 	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	map<size_t, const DMCThrown*> locIDMap; //size_t is the myid
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
		locIDMap[locMCThrowns[loc_i]->myid] = locMCThrowns[loc_i];

	locThrownSteps.clear();
	if(locMCThrowns.empty())
		return;
	locThrownSteps.push_back(pair<const DMCThrown*, deque<const DMCThrown*> >(NULL, deque<const DMCThrown*>() ) );

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(IsResonance(locMCThrowns[loc_i]->PID()))
			continue; //don't include resonances in DReaction!!

		if(locMCThrowns[loc_i]->PID() == Unknown)
			continue; //could be some weird pythia "resonance" like a diquark: just ignore them all

		//initial checks of parent id
		int locParentID = locMCThrowns[loc_i]->parentid;
		if(locParentID == 0) //photoproduced
		{
			locThrownSteps[0].second.push_back(locMCThrowns[loc_i]);
			continue;
		}
		if(locIDMap.find(locParentID) == locIDMap.end()) //produced from a particle that was not saved: spurious, don't save (e.g. product of BCAL shower)
			continue;

		//initial checks of parent pid
		Particle_t locParentPID = locIDMap[locParentID]->PID();
		bool locDoneFlag = false;
		while(((locParentPID == Unknown) || IsResonance(locParentPID)) && (!locDoneFlag))
		{
			//intermediate particle, continue towards the source
			locParentID = locIDMap[locParentID]->parentid; //parent's parent
			if(locParentID == 0) //photoproduced
			{
				locThrownSteps[0].second.push_back(locMCThrowns[loc_i]);
				locDoneFlag = true;
			}
			else if(locIDMap.find(locParentID) == locIDMap.end()) //produced from a particle that was not saved: spurious, don't save (e.g. product of BCAL shower)
				locDoneFlag = true;
			else
				locParentPID = locIDMap[locParentID]->PID();
		}

		if(Is_FinalStateParticle(locParentPID) == 1)
			continue; //e.g. parent is a final state particle (e.g. this is a neutrino from a pion decay)

		//check to see if the parent is already listed as a decaying particle //if so, add it to that step
		bool locListedAsDecayingFlag = false;
		for(size_t loc_j = 1; loc_j < locThrownSteps.size(); ++loc_j)
		{
			if(locThrownSteps[loc_j].first->myid != locParentID)
				continue;
			locThrownSteps[loc_j].second.push_back(locMCThrowns[loc_i]);
			locListedAsDecayingFlag = true;
			break;
		}
		if(locListedAsDecayingFlag)
			continue;

		//would add a new decay step, but first make sure that its parent is a decay product of a previous step
			//if the parent was not saved as a product, it may have been a decay product of a final state particle: don't save
		const DMCThrown* locThrownParent = locIDMap[locParentID];
		bool locFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locThrownSteps.size(); ++loc_j)
		{
			for(size_t loc_k = 0; loc_k < locThrownSteps[loc_j].second.size(); ++loc_k)
			{
				if(locThrownSteps[loc_j].second[loc_k] != locThrownParent)
					continue;
				locFoundFlag = true;
				break;
			}
			if(locFoundFlag)
				break;
		}
		if(!locFoundFlag)
			continue;

		//else add a new decay step and add this particle to it
		locThrownSteps.push_back(pair<const DMCThrown*, deque<const DMCThrown*> >(locThrownParent, deque<const DMCThrown*>(1, locMCThrowns[loc_i]) ));
	}

/*
cout << "THROWN STEPS: " << endl;
for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
{
	cout << ((loc_i == 0) ? 0 : locThrownSteps[loc_i].first->myid) << ": ";
	for(size_t loc_j = 0; loc_j < locThrownSteps[loc_i].second.size(); ++loc_j)
		cout << locThrownSteps[loc_i].second[loc_j]->myid << ", ";
	cout << endl;
}
*/
}

bool DAnalysisUtilities::Are_ThrownPIDsSameAsDesired(JEventLoop* locEventLoop, const deque<Particle_t>& locDesiredPIDs, Particle_t locMissingPID) const
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns, "FinalState");
	deque<Particle_t> locDesiredPIDs_Copy = locDesiredPIDs;

	bool locMissingPIDMatchedFlag = false;
	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		Particle_t locPID = (Particle_t)(locMCThrowns[loc_i]->type);

		if((!locMissingPIDMatchedFlag) && (locMissingPID == locPID))
		{
			//matched missing
			locMissingPIDMatchedFlag = true;
			continue;
		}

		bool locPIDFoundFlag = false;
		for(deque<Particle_t>::iterator locIterator = locDesiredPIDs_Copy.begin(); locIterator != locDesiredPIDs_Copy.end(); ++locIterator)
		{
			if(*locIterator != locPID)
				continue;
			locDesiredPIDs_Copy.erase(locIterator);
			locPIDFoundFlag = true;
			break;
		}
		if(!locPIDFoundFlag)
			return false;
	}

	return (locDesiredPIDs_Copy.empty());
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, unsigned int> > locSourceObjects;
	return Calc_MissingP4(locReaction, locParticleCombo, 0, -1, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	return Calc_MissingP4(locReaction, locParticleCombo, 0, -1, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, unsigned int> > locSourceObjects;
	return Calc_MissingP4(locReaction, locParticleCombo, locStepIndex, locUpToStepIndex, locUpThroughIndices, locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_MissingP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	//NOTE: this routine assumes that the p4 of a charged decaying particle with a detached vertex is the same at both vertices!
	//assumes missing particle is not the beam particle
	if(locUseKinFitDataFlag && (locParticleCombo->Get_KinFitResults() == NULL))
		return Calc_MissingP4(locReaction, locParticleCombo, locStepIndex, locUpToStepIndex, locUpThroughIndices, locSourceObjects, false); //kinematic fit failed

	DLorentzVector locMissingP4;
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	auto locReactionStep = locReaction->Get_ReactionStep(locStepIndex);

	const DKinematicData* locKinematicData = NULL;
	if(locStepIndex == 0)
	{
		//initial particle
		locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		locSourceObjects.insert(pair<const JObject*, unsigned int>(locKinematicData, abs(PDGtype(locKinematicData->PID())))); //want to use source objects for comparing
		if(locUseKinFitDataFlag) //kinfit
			locKinematicData = locParticleComboStep->Get_InitialParticle();
		locMissingP4 += locKinematicData->lorentzMomentum();
	}

	//target particle
	Particle_t locPID = locReactionStep->Get_TargetPID();
	if(locPID != Unknown)
	{
		double locMass = ParticleMass(locPID);
		locMissingP4 += DLorentzVector(DVector3(0.0, 0.0, 0.0), locMass);
	}

	auto locParticles = locUseKinFitDataFlag ? locParticleComboStep->Get_FinalParticles() : locParticleComboStep->Get_FinalParticles_Measured();
	for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
	{
		if(int(loc_j) == locReactionStep->Get_MissingParticleIndex())
			continue; //exclude missing particle
		if((int(locStepIndex) == locUpToStepIndex) && (locUpThroughIndices.find(loc_j) == locUpThroughIndices.end()))
			continue; //skip it: don't want to include it

		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_j);
		if(locDecayStepIndex > 0) //decaying-particle
		{
			//why plus? because the minus-signs are already applied during the call below
			locMissingP4 += Calc_MissingP4(locReaction, locParticleCombo, locDecayStepIndex, locUpToStepIndex, locUpThroughIndices, locSourceObjects, locUseKinFitDataFlag); //p4 returned is already < 0
		}
		else //detected
		{
			Particle_t locPID = locReactionStep->Get_FinalPID(loc_j);
			auto locDetectedP4 = locParticles[loc_j]->lorentzMomentum();
			locMissingP4 -= locDetectedP4;
			locSourceObjects.insert(pair<const JObject*, unsigned int>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_j), abs(PDGtype(locPID))));
		}
	}

	return locMissingP4;
}

TMatrixFSym DAnalysisUtilities::Calc_MissingP3Covariance(const DReaction* locReaction, const DParticleCombo* locParticleCombo) const
{
	//uses measured data!
	return Calc_MissingP3Covariance(locReaction, locParticleCombo, 0, -1, set<size_t>());
}

TMatrixFSym DAnalysisUtilities::Calc_MissingP3Covariance(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, int locUpToStepIndex, set<size_t> locUpThroughIndices) const
{
	//uses measured data!

	//NOTE: this routine assumes that the p4 of a charged decaying particle with a detached vertex is the same at both vertices!
	//assumes missing particle is not the beam particle

	//missing covariance is just sum of covariance matrices of all particles used in the calculation
		//because errors are uncorrelated: this doesn't work on kinfit data: just use kinfit matrix from missing particle then
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	auto locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	TMatrixFSym locMissingCovarianceMatrix(3);
	locMissingCovarianceMatrix.Zero();

	const DKinematicData* locKinematicData = NULL;
	if(locStepIndex == 0)
	{
		//initial particle
		locKinematicData = locParticleComboStep->Get_InitialParticle_Measured();
		TMatrixFSym locParticleCovarianceMatrix = *(locKinematicData->errorMatrix().get());
		locParticleCovarianceMatrix.ResizeTo(3, 3);
		locMissingCovarianceMatrix += locParticleCovarianceMatrix;
	}

	auto locParticles = locParticleComboStep->Get_FinalParticles_Measured();
	for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
	{
		if(int(loc_j) == locReactionStep->Get_MissingParticleIndex())
			continue; //exclude missing particle
		if((int(locStepIndex) == locUpToStepIndex) && (locUpThroughIndices.find(loc_j) == locUpThroughIndices.end()))
			continue; //skip it: don't want to include it

		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_j);
		if(locDecayStepIndex > 0) //decaying-particle
			locMissingCovarianceMatrix += Calc_MissingP3Covariance(locReaction, locParticleCombo, locDecayStepIndex, locUpToStepIndex, locUpThroughIndices);
		else //detected
		{
			TMatrixFSym locParticleCovarianceMatrix = *(locParticles[loc_j]->errorMatrix().get());
			locParticleCovarianceMatrix.ResizeTo(3, 3);
			locMissingCovarianceMatrix += locParticleCovarianceMatrix;
		}
	}

	return locMissingCovarianceMatrix;
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, unsigned int> > locSourceObjects;
	return Calc_FinalStateP4(locReaction, locParticleCombo, locStepIndex, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	return Calc_FinalStateP4(locReaction, locParticleCombo, locStepIndex, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<size_t> locToIncludeIndices, bool locUseKinFitDataFlag) const
{
	set<pair<const JObject*, unsigned int> > locSourceObjects;
	return Calc_FinalStateP4(locReaction, locParticleCombo, locStepIndex, locToIncludeIndices, locSourceObjects, locUseKinFitDataFlag);
}

DLorentzVector DAnalysisUtilities::Calc_FinalStateP4(const DReaction* locReaction, const DParticleCombo* locParticleCombo, size_t locStepIndex, set<size_t> locToIncludeIndices, set<pair<const JObject*, unsigned int> >& locSourceObjects, bool locUseKinFitDataFlag) const
{
	if(locUseKinFitDataFlag && (locParticleCombo->Get_KinFitResults() == NULL))
		return Calc_FinalStateP4(locReaction, locParticleCombo, locStepIndex, locToIncludeIndices, locSourceObjects, false); //kinematic fit failed

	DLorentzVector locFinalStateP4;
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	if(locParticleComboStep == NULL)
		return (DLorentzVector());
	auto locReactionStep = locReaction->Get_ReactionStep(locStepIndex);

	auto locParticles = locUseKinFitDataFlag ? locParticleComboStep->Get_FinalParticles() : locParticleComboStep->Get_FinalParticles_Measured();

	//subtract rescattering target if any!!
	if(locStepIndex != 0)
	{
		Particle_t locPID = locReactionStep->Get_TargetPID();
		if(locPID != Unknown)
			locFinalStateP4 -= DLorentzVector(DVector3(0.0, 0.0, 0.0), ParticleMass(locPID));
	}

	bool locDoSubsetFlag = !locToIncludeIndices.empty();
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(locDoSubsetFlag && (locToIncludeIndices.find(loc_i) == locToIncludeIndices.end()))
			continue; //skip it: don't want to include it

		if(locReactionStep->Get_MissingParticleIndex() == int(loc_i))
			return (DLorentzVector()); //bad!

		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		if(locDecayStepIndex >= 0) //decaying particle
		{
			//measured results, or not constrained by kinfit (either non-fixed mass or excluded from kinfit)
			if((!locUseKinFitDataFlag) || (!IsFixedMass(locReactionStep->Get_FinalPID(loc_i))))
				locFinalStateP4 += Calc_FinalStateP4(locReaction, locParticleCombo, locDecayStepIndex, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
			else //want kinfit results, and decaying particle p4 is constrained by kinfit
			{
				locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
				//still need source objects of decay products! dive down anyway, but ignore p4 result
				Calc_FinalStateP4(locReaction, locParticleCombo, locDecayStepIndex, set<size_t>(), locSourceObjects, locUseKinFitDataFlag);
			}
		}
		else //detected particle
		{
			locFinalStateP4 += locParticles[loc_i]->lorentzMomentum();
			locSourceObjects.insert(pair<const JObject*, unsigned int>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i), abs(PDGtype(locParticles[loc_i]->PID()))));
		}
	}
	return locFinalStateP4;
}

double DAnalysisUtilities::Calc_Energy_UnusedShowers(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo) const
{
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);
	const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();
	double locRFTime = (locEventRFBunch != NULL) ? locEventRFBunch->dTime : numeric_limits<double>::quiet_NaN();

	vector<const DNeutralShower*> locUnusedNeutralShowers;
	Get_UnusedNeutralShowers(locEventLoop, locParticleCombo, locUnusedNeutralShowers);
	
	double locEnergy_UnusedShowers = 0.; 
	for(size_t loc_i = 0; loc_i < locUnusedNeutralShowers.size(); ++loc_i) {
		const DNeutralShower* locUnusedNeutralShower = locUnusedNeutralShowers[loc_i];

		// requirements on unused showers
		double locFlightTime = (locUnusedNeutralShower->dSpacetimeVertex.Vect() - locVertex).Mag()/SPEED_OF_LIGHT;
		double locDeltaT = locUnusedNeutralShower->dSpacetimeVertex.T() - locFlightTime - locRFTime;
		double locDetectorTheta = (locUnusedNeutralShower->dSpacetimeVertex.Vect()-locVertex).Theta()*180./TMath::Pi();
		if(locDetectorTheta < 2.0 || fabs(locDeltaT) > 4.)
			continue;		

		locEnergy_UnusedShowers += locUnusedNeutralShower->dEnergy;
	}
	
	return locEnergy_UnusedShowers;
}

int DAnalysisUtilities::Calc_Momentum_UnusedTracks(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, double &locSumPMag_UnusedTracks, TVector3 &locSumP3_UnusedTracks) const
{
	vector<const DChargedTrack*> locUnusedChargedTracks;
	Get_UnusedChargedTracks(locEventLoop, locParticleCombo, locUnusedChargedTracks);
	
	for(size_t loc_i = 0; loc_i < locUnusedChargedTracks.size(); ++loc_i) {
		const DChargedTrack* locUnusedChargedTrack = locUnusedChargedTracks[loc_i];
		const DChargedTrackHypothesis *locUnusedChargedTrackHypothesis = locUnusedChargedTrack->Get_BestTrackingFOM();

		locSumPMag_UnusedTracks += locUnusedChargedTrackHypothesis->pmag();
		locSumP3_UnusedTracks += locUnusedChargedTrackHypothesis->momentum();
	}
	
	return (int)locUnusedChargedTracks.size();
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinematicData, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinematicData* locKinematicData, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir = (1.0/locKinematicData->momentum().Mag())*locKinematicData->momentum();
	DVector3 locPosition = locKinematicData->position();
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex) const
{
	DVector3 locPOCA;
	return Calc_DOCAToVertex(locKinFitParticle, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DKinFitParticle* locKinFitParticle, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locUnitDir(locKinFitParticle->Get_Momentum().Unit().X(),locKinFitParticle->Get_Momentum().Unit().Y(),locKinFitParticle->Get_Momentum().Unit().Z());
	DVector3 locPosition(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z());
	return Calc_DOCAToVertex(locUnitDir, locPosition, locVertex, locPOCA);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAToVertex(const DVector3& locUnitDir, const DVector3& locPosition, const DVector3& locVertex, DVector3& locPOCA) const
{
	DVector3 locPOCA2;
	return Calc_DOCA(locUnitDir, locUnitDir, locPosition, locVertex, locPOCA, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3& locDOCAVertex) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCAVertex(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locDOCAVertex);
}

double DAnalysisUtilities::Calc_DOCAVertex(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3& locDOCAVertex) const
{
	DVector3 locPOCA1, locPOCA2;
	double locDOCA = Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
	locDOCAVertex = 0.5*(locPOCA1 + locPOCA2);
	return locDOCA;
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinFitParticle1, locKinFitParticle2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locKinematicData1, locKinematicData2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinFitParticle* locKinFitParticle1, const DKinFitParticle* locKinFitParticle2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1(locKinFitParticle1->Get_Momentum().Unit().X(),locKinFitParticle1->Get_Momentum().Unit().Y(),locKinFitParticle1->Get_Momentum().Unit().Z());
	DVector3 locUnitDir2(locKinFitParticle2->Get_Momentum().Unit().X(),locKinFitParticle2->Get_Momentum().Unit().Y(),locKinFitParticle2->Get_Momentum().Unit().Z());
	DVector3 locVertex1(locKinFitParticle1->Get_Position().X(),locKinFitParticle1->Get_Position().Y(),locKinFitParticle1->Get_Position().Z());
	DVector3 locVertex2(locKinFitParticle2->Get_Position().X(),locKinFitParticle2->Get_Position().Y(),locKinFitParticle2->Get_Position().Z());
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DKinematicData* locKinematicData1, const DKinematicData* locKinematicData2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	DVector3 locUnitDir1 = (1.0/locKinematicData1->momentum().Mag())*locKinematicData1->momentum();
	DVector3 locUnitDir2 = (1.0/locKinematicData2->momentum().Mag())*locKinematicData2->momentum();
	DVector3 locVertex1 = locKinematicData1->position();
	DVector3 locVertex2 = locKinematicData2->position();
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2) const
{
	DVector3 locPOCA1, locPOCA2;
	return Calc_DOCA(locUnitDir1, locUnitDir2, locVertex1, locVertex2, locPOCA1, locPOCA2);
}

double DAnalysisUtilities::Calc_DOCA(const DVector3 &locUnitDir1, const DVector3 &locUnitDir2, const DVector3 &locVertex1, const DVector3 &locVertex2, DVector3 &locPOCA1, DVector3 &locPOCA2) const
{
	//originated from code by Jorn Langheinrich
	//you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
	double locUnitDot = locUnitDir1.Dot(locUnitDir2);
	double locDenominator = locUnitDot*locUnitDot - 1.0; // scalar product of directions
	double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point

	if(fabs(locDenominator) < 1.0e-15) //parallel
	{
		locDistVertToInterDOCA1 = (locVertex2 - locVertex1).Dot(locUnitDir2)/locUnitDot; //the opposite
		locDistVertToInterDOCA2 = (locVertex1 - locVertex2).Dot(locUnitDir1)/locUnitDot;
	}
	else
	{
		double locA = (locVertex1 - locVertex2).Dot(locUnitDir1);
		double locB = (locVertex1 - locVertex2).Dot(locUnitDir2);
		locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
		locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
	}

	locPOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
	locPOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
	return (locPOCA1 - locPOCA2).Mag();
}

double DAnalysisUtilities::Calc_CrudeTime(const vector<const DKinematicData*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	DVector3 locMomentum;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
		locDeltaVertex = locPOCA - locParticles[loc_i]->position();
		locMomentum = locParticles[loc_i]->momentum();
		double locTime = locParticles[loc_i]->time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->energy()/(29.9792458*locMomentum.Mag2());
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

double DAnalysisUtilities::Calc_CrudeTime(const vector<DKinFitParticle*>& locParticles, const DVector3& locCommonVertex) const
{
	//crudely propagate the track times to the common vertex and return the average track time
	DVector3 locPOCA;
	DVector3 locDeltaVertex;
	double locAverageTime = 0.0;
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		double locTime = 0.0;
		double locE = locParticles[loc_i]->Get_ShowerEnergy();
		if((locParticles[loc_i]->Get_Charge() == 0) && (locE > 0.0))
		{
			double locMass = locParticles[loc_i]->Get_Mass();
			double locPMag = sqrt(locE*locE - locMass*locMass);
			DVector3 locPosition = locParticles[loc_i]->Get_Position();
			DVector3 locDPosition(locPosition.X(), locPosition.Y(), locPosition.Z());
			DVector3 locDeltaVertex = locDPosition - locCommonVertex;
			locTime = locParticles[loc_i]->Get_Time() + locDeltaVertex.Mag()*locE/(29.9792458*locPMag);
		}
		else
		{
			Calc_DOCAToVertex(locParticles[loc_i], locCommonVertex, locPOCA);
			locDeltaVertex = locPOCA - DVector3(locParticles[loc_i]->Get_Position().X(),locParticles[loc_i]->Get_Position().Y(),locParticles[loc_i]->Get_Position().Z());
			DVector3 locMomentum(locParticles[loc_i]->Get_Momentum().X(),locParticles[loc_i]->Get_Momentum().Y(),locParticles[loc_i]->Get_Momentum().Z());
			locTime = locParticles[loc_i]->Get_Time() + locDeltaVertex.Dot(locMomentum)*locParticles[loc_i]->Get_Energy()/(29.9792458*locMomentum.Mag2());
		}
		locAverageTime += locTime;
	}
	return locAverageTime/(double(locParticles.size()));
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const vector<const DTrackTimeBased*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const vector<const DChargedTrackHypothesis*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const vector<const DKinematicData*>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;
	if(locParticles.size() == 1)
		return locParticles[0]->position();

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j], locParticles[loc_k], locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

DVector3 DAnalysisUtilities::Calc_CrudeVertex(const vector<shared_ptr<DKinFitParticle>>& locParticles) const
{
	//assumes tracks are straight lines
	//uses the midpoint of the smallest DOCA line
	DVector3 locVertex(0.0, 0.0, dTargetZCenter);

	if(locParticles.size() == 0)
		return locVertex;

	if(locParticles.size() == 1)
	  return DVector3(locParticles[0]->Get_Position().X(), locParticles[0]->Get_Position().Y(), locParticles[0]->Get_Position().Z());

	double locDOCA, locSmallestDOCA;
	DVector3 locTempVertex;

	locSmallestDOCA = 9.9E9;
	for(int loc_j = 0; loc_j < (int(locParticles.size()) - 1); ++loc_j)
	{
		for(size_t loc_k = loc_j + 1; loc_k < locParticles.size(); ++loc_k)
		{
			locDOCA = Calc_DOCAVertex(locParticles[loc_j].get(), locParticles[loc_k].get(), locTempVertex);
			if(locDOCA < locSmallestDOCA)
			{
				locSmallestDOCA = locDOCA;
				locVertex = locTempVertex;
			}
		}
	}
	return locVertex;
}

set<set<size_t> > DAnalysisUtilities::Build_IndexCombos(const DReactionStep* locReactionStep, deque<Particle_t> locToIncludePIDs) const
{
	//if locToIncludePIDs is empty, will return one set with all (except missing)
	set<set<size_t> > locCombos;

	auto locReactionStepPIDs = locReactionStep->Get_FinalPIDs();
	int locMissingParticleIndex = locReactionStep->Get_MissingParticleIndex();

	if(locToIncludePIDs.empty())
	{
		set<size_t> locCombo;
		for(size_t loc_i = 0; loc_i < locReactionStepPIDs.size(); ++loc_i)
		{
			if(int(loc_i) == locMissingParticleIndex)
				continue;
			locCombo.insert(loc_i);
		}
		locCombos.insert(locCombo);
		return locCombos;
	}

	//deque indices corresponds to locToIncludePIDs, and each set is what could pass for it
	deque<deque<size_t> > locPossibilities(locToIncludePIDs.size(), deque<size_t>());
	deque<int> locResumeAtIndices(locToIncludePIDs.size(), 0);

	//build possibilities: loop over reaction PIDs
	for(size_t loc_i = 0; loc_i < locReactionStepPIDs.size(); ++loc_i)
	{
		if(int(loc_i) == locMissingParticleIndex)
			continue;
		Particle_t locPID = locReactionStepPIDs[loc_i];

		//register where this is a valid option: loop over to-include PIDs
		for(size_t loc_j = 0; loc_j < locToIncludePIDs.size(); ++loc_j)
		{
			if(locToIncludePIDs[loc_j] == locPID)
				locPossibilities[loc_j].push_back(loc_i);
		}
	}

	//build combos
	int locParticleIndex = 0;
	deque<size_t> locComboDeque;
	while(true)
	{
		if(locParticleIndex == int(locPossibilities.size())) //end of combo: save it
		{
			set<size_t> locComboSet; //convert set to deque
			for(size_t loc_i = 0; loc_i < locComboDeque.size(); ++loc_i)
				locComboSet.insert(locComboDeque[loc_i]);
			locCombos.insert(locComboSet); //saved

			if(!Handle_Decursion(locParticleIndex, locComboDeque, locResumeAtIndices, locPossibilities))
				break;
			continue;
		}

		int& locResumeAtIndex = locResumeAtIndices[locParticleIndex];

		//if two identical pids: locResumeAtIndex must always be >= the previous locResumeAtIndex (prevents duplicates) e.g. g, p -> p, pi0, pi0
		//search for same pid previously in this step
		Particle_t locToIncludePID = locToIncludePIDs[locParticleIndex];
		for(int loc_i = locParticleIndex - 1; loc_i >= 0; --loc_i)
		{
			if(locToIncludePIDs[loc_i] == locToIncludePID)
			{
				if(locResumeAtIndex < locResumeAtIndices[loc_i])
					locResumeAtIndex = locResumeAtIndices[loc_i];
				break; //dupe type in step: resume-at advanced to next
			}
		}

		if(locResumeAtIndex >= int(locPossibilities[locParticleIndex].size()))
		{
			if(!Handle_Decursion(locParticleIndex, locComboDeque, locResumeAtIndices, locPossibilities))
				break;
			continue;
		}

		// pid found
		locComboDeque.push_back(locPossibilities[locParticleIndex][locResumeAtIndex]);
		++locResumeAtIndex;
		++locParticleIndex;
	}

	return locCombos;
}

bool DAnalysisUtilities::Handle_Decursion(int& locParticleIndex, deque<size_t>& locComboDeque, deque<int>& locResumeAtIndices, deque<deque<size_t> >& locPossibilities) const
{
	do
	{
		if(locParticleIndex < int(locResumeAtIndices.size())) //else just saved a combo
			locResumeAtIndices[locParticleIndex] = 0; //finding this particle failed: reset

		--locParticleIndex; //go to previous particle
		if(locParticleIndex < 0)
			return false; //end of particles: end of finding all combos

		locComboDeque.pop_back(); //reset this index
	}
	while(locResumeAtIndices[locParticleIndex] == int(locPossibilities[locParticleIndex].size()));

	return true;
}

//The POCA cannot technically be solved analytically, but we can approximate it pretty accurately
	//First, propagate the track to somewhere very close to the true POCA
	//Then, the equation can be solved by substituting cos(x) = 1, sin(x) = x
double DAnalysisUtilities::Propagate_Track(int locCharge, const DVector3& locPropagateToPoint, DLorentzVector& locMeasuredX4, DLorentzVector& locMeasuredP4, TMatrixFSym* locCovarianceMatrix) const
{
	//ASSUMES THAT THE B-FIELD IS IN THE +Z DIRECTION!!!!!!!!!!!!!!!
	double locDistance = (locMeasuredX4.Vect() - locPropagateToPoint).Mag();
	if(!Get_IsBFieldNearBeamline() || (locCharge == 0) || (locDistance < dMinDistanceForStraightTrack))
	{
		//use simpler methods
		DVector3 locPOCA;
		auto locP3 = locMeasuredP4.Vect();
		Calc_DOCAToVertex(locP3.Unit(), locMeasuredX4.Vect(), locPropagateToPoint, locPOCA);
		locMeasuredX4.SetVect(locPOCA);
		auto locDistanceVector = locPOCA - locMeasuredX4.Vect();
		//negative: if you had to propagate the track forwards, that means the path length is LESS than what you thought it was
		auto locDeltaPathLength = (locP3.Dot(locDistanceVector) > 0.0) ? -1.0*locDistanceVector.Mag() : locDistanceVector.Mag();
		locMeasuredX4.SetT(locMeasuredX4.T() - locDeltaPathLength/(locMeasuredP4.Beta()*SPEED_OF_LIGHT)); //v = s/t, t = s/v = s/(beta*c) = s*E/(p*c) =
		return locDeltaPathLength;
	}

	//propagate the track to the same z as the vertex (if pz != 0)
	double locTotalDeltaPathLength = 0.0;
	DLorentzVector locTempPropagatedPosition = locMeasuredX4;
	DLorentzVector locTempPropagatedMomentum = locMeasuredP4;
	if(fabs(locMeasuredP4.Pz()) > 0.0)
	{
		locTotalDeltaPathLength += (locPropagateToPoint.Z() - locMeasuredX4.Z())*locMeasuredP4.P()/locMeasuredP4.Pz(); //s = (z - z0)*p/p0z
		Propagate_Track(locTotalDeltaPathLength, locCharge, locTempPropagatedPosition, locTempPropagatedMomentum, NULL);
	}

	//now, step along the track until we get close to the POCA
	locTotalDeltaPathLength += Calc_PathLength_Step(locCharge, locPropagateToPoint, locTempPropagatedPosition, locTempPropagatedMomentum);

	//now, the path length is very close to accurate. get the rest of the way there
	locTotalDeltaPathLength += Calc_PathLength_FineGrained(locCharge, locPropagateToPoint, locTempPropagatedPosition.Vect(), locTempPropagatedMomentum.Vect());

	//Finally, propagate the track the given distance and return
//cout << "FINAL: xyz, distance = " << locMeasuredX4.X() << ", " << locMeasuredX4.Y() << ", " << locMeasuredX4.Z() << ", " << locTotalDeltaPathLength << endl;
	Propagate_Track(locTotalDeltaPathLength, locCharge, locMeasuredX4, locMeasuredP4, locCovarianceMatrix);

	//negative: if you had to propagate the track forwards, that means the path length is LESS than what you thought it was
	return -1.0*locTotalDeltaPathLength;
}

double DAnalysisUtilities::Calc_PathLength_Step(int locCharge, const DVector3& locPropagateToPoint, DLorentzVector& locMeasuredX4, DLorentzVector& locMeasuredP4) const
{
	//now, step slowly along the track (in path length), trying to get closer to the POCA
		//for the final calculation to work, rho*s must be small: cos(rho*s) -> 1
	DVector3 locBField = Get_BField(locPropagateToPoint);
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	double locPMag = locMeasuredP4.P();
	double locPathLengthOneRotation = 2.0*TMath::Pi()*locPMag/locA;
	double locStepDistance = locPathLengthOneRotation/12.0; //30 degrees //e.g. for 90 degree tracks

	//guard against long steps in z (e.g. for 2 degree tracks)
	double locStepDeltaZ = locMeasuredP4.Pz()*locStepDistance/locPMag;
	double locTargetDeltaZ = 1.0;
	if(fabs(locStepDeltaZ) > 5.0)
		locStepDistance = locTargetDeltaZ*locPMag/locMeasuredP4.Pz();

	//First, try stepping in +pvec direction
	double locStepDirection = 1.0;
	double locTotalDeltaPathLength = 0.0;
	double locDeltaX3 = (locMeasuredX4.Vect() - locPropagateToPoint).Mag();
	double locRhoS = locStepDistance*locA/locPMag;
	double locRhoSCutOff = 0.1;
	while(locRhoS > locRhoSCutOff) //at this point, cos(rho*s) = 0.995: close enough
	{
		double locDeltaPathLength = locStepDirection*locStepDistance;
		locTotalDeltaPathLength += locDeltaPathLength;
//cout << "STEP LOOP: xyz, distance = " << locMeasuredX4.X() << ", " << locMeasuredX4.Y() << ", " << locMeasuredX4.Z() << ", " << locDeltaPathLength << endl;
		Propagate_Track(locDeltaPathLength, locCharge, locMeasuredX4, locMeasuredP4, NULL);

		double locNewDeltaX3 = (locMeasuredX4.Vect() - locPropagateToPoint).Mag();
		bool locGettingCloserFlag = (locNewDeltaX3 < locDeltaX3);
		locDeltaX3 = locNewDeltaX3;

		if(locGettingCloserFlag)
			continue; //stepping in correct direction, can still try to get closer

		//are now farther away than before: reverse direction, decrease step size
		locStepDirection *= -1.0;
		locStepDistance /= 2.0;
		locRhoS = locStepDistance*locA/locPMag;

		//if about to break, revert to older, better position
		if(locRhoS < locRhoSCutOff)
		{
			locTotalDeltaPathLength -= locDeltaPathLength;
			Propagate_Track(-1.0*locDeltaPathLength, locCharge, locMeasuredX4, locMeasuredP4, NULL);
		}
	}

	return locTotalDeltaPathLength;
}

double DAnalysisUtilities::Calc_PathLength_FineGrained(int locCharge, const DVector3& locPropagateToPoint, DVector3 locMeasuredPosition, DVector3 locMeasuredMomentum) const
{
	//ASSUMES B-FIELD IS IN +Z DIRECTION!!

	//distance = sqrt((delta-x)^2 + (delta-y)^2 + (delta-z)^2)
	//find the delta-path-length s that minimizes the distance equation
	//take derivative of the distance equation, set = 0

	//Mathematica yields: F*(B*s + C) + D*cos(rho*s) + E*sin(rho*s) = 0
	//(where BCDEF are various, non-s-dependent terms)
	//This is unsolvable analytically. However, we know that s is small, since we're almost there

	//So, substitute cos(x) = 1, sin(x) = x:
	//F*(B*s + C) + D + E*rho*s = 0
	//Solving gives: s = - (F*C + D) / (F*B + E*rho)
	//Note that B = F

	DVector3 locBField = Get_BField(locMeasuredPosition);
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
	double locPMag = locMeasuredMomentum.Mag();
	double locRho = locA/locPMag;

	DVector3 locDeltaX3 = locMeasuredPosition - locPropagateToPoint;
	double locF = locMeasuredMomentum.Pz()/locPMag;
	double locC = locDeltaX3.Z();
	double locD = locRho*(locMeasuredMomentum.Px()*locDeltaX3.X() + locMeasuredMomentum.Py()*locDeltaX3.Y())/locA;
	double locPerpPSq = locMeasuredMomentum.Px()*locMeasuredMomentum.Px() + locMeasuredMomentum.Py()*locMeasuredMomentum.Py();
	double locE = locRho*(locPerpPSq/locA - locMeasuredMomentum.Py()*locDeltaX3.X() + locMeasuredMomentum.Px()*locDeltaX3.Y())/locA;
	return -1.0*(locF*locC + locD)/(locF*locF + locE*locRho);
}

void DAnalysisUtilities::Propagate_Track(double locDeltaPathLength, int locCharge, DLorentzVector& locX4, DLorentzVector& locP4, TMatrixFSym* locCovarianceMatrix) const
{
	//ASSUMES THAT THE B-FIELD IS IN THE +Z DIRECTION!!!!!!!!!!!!!!!

	DVector3 locBField = Get_BField(locX4.Vect());
	if(!(locBField.Mag() > 0.0))
		return;
	DVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	double locPMag = locP4.P();
	double locRhoS = locDeltaPathLength*locA/locPMag;

	DLorentzVector locDeltaX4; //x - x0
	locDeltaX4.SetX(locP4.Px()*sin(locRhoS)/locA - locP4.Py()*(1.0 - cos(locRhoS))/locA);
	locDeltaX4.SetY(locP4.Py()*sin(locRhoS)/locA + locP4.Px()*(1.0 - cos(locRhoS))/locA);
	locDeltaX4.SetZ(locP4.Pz()*locDeltaPathLength/locPMag);
	locDeltaX4.SetT(locDeltaPathLength/(locP4.Beta()*SPEED_OF_LIGHT)); //v = s/t, t = s/v = s/(beta*c) = s*E/(p*c) =
	locX4 += locDeltaX4;

	if(locCovarianceMatrix == NULL)
	{
		locP4.SetVect(locP4.Vect() - locDeltaX4.Vect().Cross(locA*locH));
		return; //don't update error matrix
	}

	//transform(i, j) = di/dj, i = new, j = old //pxyz, xyz, t
	TMatrixF locTransformMatrix(7, 7);
	locTransformMatrix.Zero();

	double locPCubed = locPMag*locPMag*locPMag;
	double locSOverP3 = locDeltaPathLength/locPCubed;

	double locSOverP3_PxPx = locSOverP3*locP4.Px()*locP4.Px();
	double locSOverP3_PxPy = locSOverP3*locP4.Px()*locP4.Py();
	double locSOverP3_PxPz = locSOverP3*locP4.Px()*locP4.Pz();

	double locSOverP3_PyPy = locSOverP3*locP4.Py()*locP4.Py();
	double locSOverP3_PyPz = locSOverP3*locP4.Py()*locP4.Pz();
	double locSOverP3_PzPz = locSOverP3*locP4.Pz()*locP4.Pz();

	//dpx
	locTransformMatrix(0, 0) = (1.0 + locA*locSOverP3_PxPy)*cos(locRhoS) + locA*locSOverP3_PxPx*sin(locRhoS);
	locTransformMatrix(0, 1) = (-1.0 + locA*locSOverP3_PxPy)*sin(locRhoS) + locA*locSOverP3_PyPy*cos(locRhoS);
	locTransformMatrix(0, 2) = locA*locSOverP3_PyPz*cos(locRhoS) + locA*locSOverP3_PxPz*sin(locRhoS);

	//dpy
	locTransformMatrix(1, 0) = -1.0*locA*locSOverP3_PxPx*cos(locRhoS) + (1.0 + locA*locSOverP3_PxPy)*sin(locRhoS);
	locTransformMatrix(1, 1) = (1.0 - locA*locSOverP3_PxPy)*cos(locRhoS) + locA*locSOverP3_PyPy*sin(locRhoS);
	locTransformMatrix(1, 2) = locA*locSOverP3_PyPz*sin(locRhoS) - locA*locSOverP3_PxPz*cos(locRhoS);

	//dpz
	locTransformMatrix(2, 2) = 1.0;

	//dx
	locTransformMatrix(3, 0) = (1.0/locA + locSOverP3_PxPy)*sin(locRhoS) - locSOverP3_PxPx*cos(locRhoS);
	locTransformMatrix(3, 1) = -1.0/locA + (1.0/locA - locSOverP3_PxPy)*cos(locRhoS) + locSOverP3_PyPy*sin(locRhoS);
	locTransformMatrix(3, 2) = locSOverP3_PyPz*sin(locRhoS) - locSOverP3_PxPz*cos(locRhoS);
	locTransformMatrix(3, 3) = 1.0;

	//dy
	locTransformMatrix(4, 0) = 1.0/locA - (1.0/locA + locSOverP3_PxPy)*cos(locRhoS) + locSOverP3_PxPx*sin(locRhoS);
	locTransformMatrix(4, 1) = (1.0/locA - locSOverP3_PxPy)*sin(locRhoS) - locSOverP3_PyPy*cos(locRhoS);
	locTransformMatrix(4, 2) = -1.0*locSOverP3_PyPz*cos(locRhoS) - locSOverP3_PxPz*sin(locRhoS);
	locTransformMatrix(4, 4) = 1.0;

	//dz
	locTransformMatrix(5, 0) = -1.0*locSOverP3_PxPz;
	locTransformMatrix(5, 1) = -1.0*locSOverP3_PyPz;
	locTransformMatrix(5, 2) = locDeltaPathLength/locPMag - locSOverP3_PzPz;
	locTransformMatrix(5, 5) = 1.0;

	//dt
	locTransformMatrix(6, 0) = -1.0*locP4.M2()*locSOverP3*locP4.Px()/(SPEED_OF_LIGHT*locP4.E());
	locTransformMatrix(6, 1) = -1.0*locP4.M2()*locSOverP3*locP4.Py()/(SPEED_OF_LIGHT*locP4.E());
	locTransformMatrix(6, 2) = -1.0*locP4.M2()*locSOverP3*locP4.Pz()/(SPEED_OF_LIGHT*locP4.E());
	locTransformMatrix(6, 6) = 1.0;

	//transform!!
	locCovarianceMatrix->Similarity(locTransformMatrix);

	//update p3 //must do at the end!!
	locP4.SetVect(locP4.Vect() - locDeltaX4.Vect().Cross(locA*locH));
}
