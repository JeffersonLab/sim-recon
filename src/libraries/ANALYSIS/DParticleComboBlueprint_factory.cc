#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DParticleComboBlueprint_factory.h"

//------------------
// init
//------------------
jerror_t DParticleComboBlueprint_factory::init(void)
{
	MAX_DParticleComboBlueprintStepPoolSize = 3000;

	dDebugLevel = 0;
	dMinProtonMomentum = pair<bool, double>(false, -1.0);
	dMinTrackingFOM = pair<bool, double>(false, -1.0);
	dReactionShowerSelectionTag = pair<bool, string>(false, "");
	dReactionTrackSelectionTag = pair<bool, string>(false, "");
	dHasDetectorMatchFlag = pair<bool, bool>(false, false);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleComboBlueprint_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	//BE CAREFUL: DON'T DO ANYTHING THAT REQUIRES THE brun() METHOD OF THIS FACTORY TO BE CALLED!!!!
	dTrackTimeBasedFactory_Combo = dynamic_cast<DTrackTimeBased_factory_Combo*>(locEventLoop->GetFactory("DTrackTimeBased", "Combo"));

	gPARMS->SetDefaultParameter("COMBOBLUEPRINTS:DEBUG_LEVEL", dDebugLevel);

	// In the following try-catch blocks, gPARMS->GetParameter will throw an
	// exception if the parameter doesn't exist leaving both the X.second and
	// X.first elements of the relevant variable untouched. If the parameter
	// does exist, the value is copied into X.second and the X.first value
	// gets set to true on the subsequent line.

	if(gPARMS->Exists("COMBO:MIN_PROTON_MOMENTUM"))
	{
		dMinProtonMomentum.first = true;
		gPARMS->GetParameter("COMBO:MIN_PROTON_MOMENTUM", dMinProtonMomentum.second);
	}

	if(gPARMS->Exists("COMBO:MIN_TRACKING_FOM"))
	{
		dMinTrackingFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_TRACKING_FOM", dMinTrackingFOM.second);
	}

	if(gPARMS->Exists("COMBO:REACTION_TRACK_SELECT_TAG"))
	{
		dReactionTrackSelectionTag.first = true;
		gPARMS->GetParameter("COMBO:REACTION_TRACK_SELECT_TAG", dReactionTrackSelectionTag.second);
	}

	if(gPARMS->Exists("COMBO:REACTION_SHOWER_SELECT_TAG"))
	{
		dReactionShowerSelectionTag.first = true;
		gPARMS->GetParameter("COMBO:REACTION_SHOWER_SELECT_TAG", dReactionShowerSelectionTag.second);
	}

	if(gPARMS->Exists("COMBO:HAS_DETECTOR_MATCH_FLAG"))
	{
		dHasDetectorMatchFlag.first = true;
		gPARMS->GetParameter("COMBO:HAS_DETECTOR_MATCH_FLAG", dHasDetectorMatchFlag.second);
	}

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);
	MAX_DParticleComboBlueprintStepPoolSize = 3000*locReactions.size();

	return NOERROR;
}

void DParticleComboBlueprint_factory::Get_Reactions(JEventLoop *locEventLoop, vector<const DReaction*>& locReactions) const
{
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	locReactions.clear();
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
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
}

//------------------
// evnt
//------------------
jerror_t DParticleComboBlueprint_factory::evnt(JEventLoop *locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DParticleComboBlueprint_factory::evnt()");
#endif

	Reset_Pools();
	dBlueprintStepMap.clear();
	dSavedBlueprintSteps.clear();

	vector<const DReaction*> locReactions;
	Get_Reactions(locEventLoop, locReactions);

	locEventLoop->GetSingle(dVertex);
	locEventLoop->GetSingle(dDetectorMatches);

	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		Build_ParticleComboBlueprints(locEventLoop, locReactions[loc_i]);

cout << "Event, # blues = " << locEventLoop->GetJEvent().GetEventNumber() << ", " << _data.size() << endl;
	return NOERROR;
}

jerror_t DParticleComboBlueprint_factory::Build_ParticleComboBlueprints(JEventLoop* locEventLoop, const DReaction* locReaction)
{
	string locReactionTrackSelectionTag = dReactionTrackSelectionTag.first ? dReactionTrackSelectionTag.second : locReaction->Get_ChargedTrackFactoryTag();
	string locReactionShowerSelectionTag = dReactionShowerSelectionTag.first ? dReactionShowerSelectionTag.second : locReaction->Get_NeutralShowerFactoryTag();

	vector<const DChargedTrack*> locChargedTrackVector;
	locEventLoop->Get(locChargedTrackVector, locReactionTrackSelectionTag.c_str());
	vector<const DNeutralShower*> locNeutralShowerVector;
	locEventLoop->Get(locNeutralShowerVector, locReactionShowerSelectionTag.c_str());

	deque<const JObject*> locNeutralShowerDeque;
	for(size_t loc_i = 0; loc_i < locNeutralShowerVector.size(); ++loc_i)
		locNeutralShowerDeque.push_back(static_cast<const JObject*>(locNeutralShowerVector[loc_i]));

	if(dDebugLevel > 0)
		cout << "Reaction Name, # Reaction steps = " << locReaction->Get_ReactionName() << ", " << locReaction->Get_NumReactionSteps() << endl;
	if(locReaction->Get_NumReactionSteps() == 0)
		return RESOURCE_UNAVAILABLE;

	//make sure not more than one missing particle
	size_t locNumMissingParticles = 0;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locReaction->Get_ReactionStep(loc_i)->Get_MissingParticleIndex() != -1)
			++locNumMissingParticles;
	}
	if(locNumMissingParticles > 1)
	{
		cout << "ERROR: Too many missing particles in DReaction.  No DParticleComboBlueprints generated." << endl;
		return RESOURCE_UNAVAILABLE;
	}

	//sort charged particles into +/-
	deque<const JObject*> locChargedTrackDeque_Positive;
	deque<const JObject*> locChargedTrackDeque_Negative;
	const DChargedTrack* locChargedTrack;
	//Note that a DChargedTrack object can sometimes contain both positively and negatively charged hypotheses simultaneously: sometimes the tracking flips the sign of the track
	for(size_t loc_i = 0; loc_i < locChargedTrackVector.size(); ++loc_i)
	{
		locChargedTrack = locChargedTrackVector[loc_i];
		bool locCouldBePositiveFlag = false;
		bool locCouldBeNegativeFlag = false;
		for(size_t loc_j = 0; loc_j < locChargedTrack->dChargedTrackHypotheses.size(); ++loc_j)
		{
			if(locChargedTrack->dChargedTrackHypotheses[loc_j]->charge() > 0.0)
				locCouldBePositiveFlag = true;
			else
				locCouldBeNegativeFlag = true;
		}
		if(locCouldBePositiveFlag)
			locChargedTrackDeque_Positive.push_back(static_cast<const JObject*>(locChargedTrack));
		if(locCouldBeNegativeFlag)
			locChargedTrackDeque_Negative.push_back(static_cast<const JObject*>(locChargedTrack));
	}

	if(dDebugLevel > 0)
		cout << "#+, #-, #0 particles = " << locChargedTrackDeque_Positive.size() << ", " << locChargedTrackDeque_Negative.size() << ", " << locNeutralShowerDeque.size() << endl;

	//set up combo loop
	deque<deque<int> > locResumeAtIndexDeque; //1st index is step, 2nd is particle (the initial particle, then the final particles)
	deque<deque<int> > locNumPossibilitiesDeque; //1st index is step, 2nd is particle (the initial particle, then the final particles)
	map<int, int> locInitialParticleStepFromIndexMap; //ints are: step index, production-step index
	map<pair<int, int>, int> locFinalStateDecayStepIndexMap; //ints are: step index, particle index, decay step index
	if(!Setup_ComboLoop(locReaction, locNeutralShowerDeque.size(), locChargedTrackVector.size(), locChargedTrackDeque_Positive.size(), locChargedTrackDeque_Negative.size(), locResumeAtIndexDeque, locNumPossibilitiesDeque, locInitialParticleStepFromIndexMap, locFinalStateDecayStepIndexMap))
	{
		if(dDebugLevel > 0)
			cout << "not enough detected particles with the correct charges for the event: no combos found." << endl;
		return NOERROR;
	}

	if(dDebugLevel > 10)
	{
		cout << "locResumeAtIndexDeque: ";
		for(size_t loc_i = 0; loc_i < locResumeAtIndexDeque.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locResumeAtIndexDeque[loc_i].size(); ++loc_j)
				cout << locResumeAtIndexDeque[loc_i][loc_j] << ", ";
			cout << ";";
		}
		cout << endl;

		cout << "locNumPossibilitiesDeque: ";
		for(size_t loc_i = 0; loc_i < locNumPossibilitiesDeque.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locNumPossibilitiesDeque[loc_i].size(); ++loc_j)
				cout << locNumPossibilitiesDeque[loc_i][loc_j] << ", ";
			cout << ";";
		}
		cout << endl;

		cout << "locInitialParticleStepFromIndexMap: ";
		map<int, int>::iterator locIterator = locInitialParticleStepFromIndexMap.begin();
		for(; locIterator != locInitialParticleStepFromIndexMap.end(); ++locIterator)
			cout << locIterator->first << ", " << locIterator->second << endl;

		cout << "locFinalStateDecayStepIndexMap: ";
		map<pair<int, int>, int>::iterator locPairIterator = locFinalStateDecayStepIndexMap.begin();
		for(; locPairIterator != locFinalStateDecayStepIndexMap.end(); ++locPairIterator)
			cout << locPairIterator->first.first << ", " << locPairIterator->first.second << ": " << locPairIterator->second << endl;
	}

	//find the combos!!
	dCurrentComboSourceObjects.clear();
	Find_Combos(locReaction, locNeutralShowerDeque, locChargedTrackDeque_Positive, locChargedTrackDeque_Negative, locResumeAtIndexDeque, locNumPossibilitiesDeque, locInitialParticleStepFromIndexMap, locFinalStateDecayStepIndexMap);

	if(dDebugLevel > 10)
	{
		cout << "print pointers: " << endl;
		for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
		{
			cout << "COMBO " << loc_i << endl;
			for(size_t loc_j = 0; loc_j < _data[loc_i]->Get_NumParticleComboBlueprintSteps(); ++loc_j)
			{
				cout << "Step " << loc_j << " pointers: ";
				for(size_t loc_k = 0; loc_k < _data[loc_i]->Get_ParticleComboBlueprintStep(loc_j)->Get_NumFinalParticleSourceObjects(); ++loc_k)
					cout << _data[loc_i]->Get_ParticleComboBlueprintStep(loc_j)->Get_FinalParticle_SourceObject(loc_k) << ", ";
				cout << endl;
			}
		}
	}

	return NOERROR;
}

bool DParticleComboBlueprint_factory::Setup_ComboLoop(const DReaction* locReaction, int locNumDetectedNeutralParticles, int locNumDetectedChargedParticles, int locNumDetectedPositiveParticles, int locNumDetectedNegativeParticles, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap)
{
	//setup locResumeAtIndexDeque, & locNumPossibilitiesDeque
	Particle_t locAnalysisPID;
	int locMissingParticleIndex;
	unsigned int locNumSteps = locReaction->Get_NumReactionSteps();
	locResumeAtIndexDeque.clear();
	locResumeAtIndexDeque.clear();
	int locCharge;
	int locNumNeededChargedParticles = 0, locNumNeededPositiveParticles = 0, locNumNeededNegativeParticles = 0, locNumNeededNeutralParticles = 0;

	locInitialParticleStepFromIndexMap[0] = -1;
	for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		size_t locNumFinalParticles = locReactionStep->Get_NumFinalParticleIDs();

		//setup final state particles num possibilities & resume-at index deque
		deque<int> locTempDeque(locNumFinalParticles, 0);
		locResumeAtIndexDeque.push_back(locTempDeque);
		locMissingParticleIndex = locReactionStep->Get_MissingParticleIndex();
		for(size_t loc_j = 0; loc_j < locNumFinalParticles; ++loc_j)
		{
			locAnalysisPID = locReactionStep->Get_FinalParticleID(loc_j);

			if(locMissingParticleIndex == int(loc_j))
			{
				locTempDeque[loc_j] = 1;
				continue;
			}

			// check to see if this particle has a decay that is represented in a future step
			// e.g. on Lambda in g, p -> K+, Lambda; where a later step is Lambda -> p, pi-
			int locResumeAtIndex = 0;
			int locDecayStepIndex = Grab_DecayingParticle(locAnalysisPID, locResumeAtIndex, locReaction, loc_i, loc_j);
			if(locDecayStepIndex >= 0)
			{
				if(dDebugLevel > 10)
					cout << "decaying particle" << endl;
				locTempDeque[loc_j] = 1;
				locFinalStateDecayStepIndexMap[pair<int, int>(loc_i, loc_j)] = locDecayStepIndex; //store step where this particle decays
				locInitialParticleStepFromIndexMap[locDecayStepIndex] = loc_i; //store step where this particle is produced
				continue;
			}

			//else use detected particles
			locCharge = ParticleCharge(locAnalysisPID);
			if(locCharge > 0)
			{
				++locNumNeededPositiveParticles;
				++locNumNeededChargedParticles;
				locTempDeque[loc_j] = locNumDetectedPositiveParticles;
			}
			else if(locCharge < 0)
			{
				++locNumNeededNegativeParticles;
				++locNumNeededChargedParticles;
				locTempDeque[loc_j] = locNumDetectedNegativeParticles;
			}
			else
			{
				++locNumNeededNeutralParticles;
				locTempDeque[loc_j] = locNumDetectedNeutralParticles;
			}
		}
		locNumPossibilitiesDeque.push_back(locTempDeque);
	}

	if((locNumNeededPositiveParticles > locNumDetectedPositiveParticles) || (locNumNeededNegativeParticles > locNumDetectedNegativeParticles))
		return false; //not enough particles of a given charge for the event
	if((locNumNeededChargedParticles > locNumDetectedChargedParticles) || (locNumNeededNeutralParticles > locNumDetectedNeutralParticles))
		return false; //not enough particles of a given charge for the event //#charged can fail here if a charged track has hypotheses with different charges

	//make sure decaying particles are valid: one entry per step except the first one
	for(size_t loc_i = 0; loc_i < locNumSteps; ++loc_i)
	{
		if(locInitialParticleStepFromIndexMap.find(loc_i) == locInitialParticleStepFromIndexMap.end())
			return false;
	}

	return true;
}

void DParticleComboBlueprint_factory::Find_Combos(const DReaction* locReaction, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap, map<pair<int, int>, int>& locFinalStateDecayStepIndexMap)
{
#ifdef VTRACE
	VT_TRACER("DParticleComboBlueprint_factory::Find_Combos()");
#endif
	DParticleComboBlueprint* locParticleComboBlueprint = new DParticleComboBlueprint();
	locParticleComboBlueprint->Set_Reaction(locReaction);

	int locStepIndex = locReaction->Get_NumReactionSteps() - 1;
	int locParticleIndex = 0; //final = 0 -> (#final - 1)
	DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
	locParticleComboBlueprintStep->Set_ReactionStep(locReaction->Get_ReactionStep(locStepIndex));
	locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(locInitialParticleStepFromIndexMap[locStepIndex]);

	do
	{
		if(dDebugLevel > 10)
			cout << "do loop: step & particle indices = " << locStepIndex << ", " << locParticleIndex << endl;
		if(locParticleIndex == int(locNumPossibilitiesDeque[locStepIndex].size()))
		{
			if(dDebugLevel > 10)
				cout << "end of reaction step" << endl;
			if(!Handle_EndOfReactionStep(locReaction, locParticleComboBlueprint, locParticleComboBlueprintStep, locStepIndex, locParticleIndex, locResumeAtIndexDeque, locNumPossibilitiesDeque, locInitialParticleStepFromIndexMap))
				break;
			continue;
		}

		//get Analysis PID & the resume-at index
		int locMissingParticleIndex = locReaction->Get_ReactionStep(locStepIndex)->Get_MissingParticleIndex();
		Particle_t locAnalysisPID = locReaction->Get_ReactionStep(locStepIndex)->Get_FinalParticleID(locParticleIndex);
		int& locResumeAtIndex = locResumeAtIndexDeque[locStepIndex][locParticleIndex];

		if(dDebugLevel > 10)
			cout << "do loop: locAnalysisPID, locMissingParticleIndex, locResumeAtIndex = " << ParticleType(locAnalysisPID) << ", " << locMissingParticleIndex << ", " << locResumeAtIndex << endl;

		//handle if this is a missing particle
		if(locMissingParticleIndex == locParticleIndex)
		{
			if(dDebugLevel > 10)
				cout << "missing particle" << endl;
			// e.g. on neutron in g, p -> pi+, (n)
			//only one possibility ("missing"), so just set NULL and advance
			locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(NULL, -1); //missing
			locResumeAtIndex = 1;
			++locParticleIndex;
			continue;
		}

		// check to see if this particle has a decay that is represented in a future step
		// e.g. on Lambda in g, p -> K+, Lambda; where a later step is Lambda -> p, pi-
		pair<int, int> locParticlePair(locStepIndex, locParticleIndex);
		if(locFinalStateDecayStepIndexMap.find(locParticlePair) != locFinalStateDecayStepIndexMap.end())
		{
			int locDecayStepIndex = locFinalStateDecayStepIndexMap[locParticlePair];
			if(dDebugLevel > 10)
				cout << "decaying particle" << endl;
			locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(NULL, locDecayStepIndex); //decaying
			locResumeAtIndex = 1;
			++locParticleIndex;
			continue;
		}

		//if two detected particles of same type in a step: locResumeAtIndex must always be >= the previous locResumeAtIndex (prevents duplicates) e.g. g, d -> p, p, pi-
		//search for same pid previously in this step (and non-missing)
		for(int loc_i = locParticleIndex - 1; loc_i >= 0; --loc_i)
		{
			if(loc_i == locMissingParticleIndex)
				continue;
			if(locReaction->Get_ReactionStep(locStepIndex)->Get_FinalParticleID(loc_i) == locAnalysisPID)
			{
				if(locResumeAtIndex < locResumeAtIndexDeque[locStepIndex][loc_i])
					locResumeAtIndex = locResumeAtIndexDeque[locStepIndex][loc_i];
				if(dDebugLevel > 10)
					cout << "dupe type in step; locResumeAtIndex = " << locResumeAtIndex << endl;
				break;
			}
		}

		// else grab a detected track
		const JObject* locSourceObject = Grab_DetectedTrack(locReaction, locAnalysisPID, locResumeAtIndex, locNeutralShowerDeque, locChargedTrackDeque_Positive, locChargedTrackDeque_Negative);
		if(locSourceObject == NULL)
		{
			if(dDebugLevel > 10)
				cout << "can't find detected particle" << endl;
			if(!Handle_Decursion(locParticleComboBlueprint, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleIndex, locStepIndex, locParticleComboBlueprintStep))
				break;
			continue;
		}

		if(dDebugLevel > 10)
			cout << "detected track found, locResumeAtIndex now = " << locResumeAtIndex << endl;
		locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(locSourceObject, -2); //detected
		dCurrentComboSourceObjects.insert(locSourceObject);
		++locParticleIndex;
	}
	while(true);
	delete locParticleComboBlueprint; //delete the last, extra one
}

bool DParticleComboBlueprint_factory::Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, map<int, int>& locInitialParticleStepFromIndexMap)
{
	//end of step

	//if all of the particle types in this step are identical to all of the particle types in a previously-done step (and none of either are missing), regardless of the order they are listed:
		//the locResumeAtIndex system of the new step's particles must be >= the locResumeAtIndex system of the particles in the LAST OF the matching previously-done steps (e.g. could be 3x pi0s)
			//in other words: if 2x pi0 -> g, g; CD, AB ok and AB, CD bad; but also BC, AD ok (BC system occurs after AD system)
	if(Check_IfDuplicateStepCombo(locParticleComboBlueprint, locParticleComboBlueprintStep, locStepIndex, locResumeAtIndexDeque, locNumPossibilitiesDeque)) //make sure none of the dupe particles are missing
	{
		if(dDebugLevel > 10)
			cout << "duplicate step combo" << endl;
		if(!Handle_Decursion(locParticleComboBlueprint, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleIndex, locStepIndex, locParticleComboBlueprintStep))
			return false;
		return true;
	}

	// cut on invariant mass if desired
	Particle_t locStepInitialPID = locParticleComboBlueprintStep->Get_InitialParticleID();
	double locMinInvariantMass = 1.0, locMaxInvariantMass = -1.0;
	if(locReaction->Get_InvariantMassCut(locStepInitialPID, locMinInvariantMass, locMaxInvariantMass))
	{
		DLorentzVector locP4;
		if(Calc_FinalStateP4(locReaction->Get_NumReactionSteps(), locParticleComboBlueprint, locParticleComboBlueprintStep, -1, locP4))
		{
			double locInvariantMass = locP4.M();
			if((locInvariantMass < locMinInvariantMass) || (locInvariantMass > locMaxInvariantMass))
			{
				if(dDebugLevel > 10)
					cout << "bad invariant mass" << endl;
				if(!Handle_Decursion(locParticleComboBlueprint, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleIndex, locStepIndex, locParticleComboBlueprintStep))
					return false;
				return true;
			}
		}
	}

	//step is good: advance to next step

	//first check to see if identical to a previous saved step; if so, just save the old step and recycle the current one
	map<DParticleComboBlueprintStep, DParticleComboBlueprintStep*>::iterator locStepIterator = dBlueprintStepMap.find(*locParticleComboBlueprintStep);
	if(locStepIterator != dBlueprintStepMap.end())
	{
		//identical step found, recycle current one
		Recycle_ParticleComboBlueprintStep(locParticleComboBlueprintStep);
		locParticleComboBlueprintStep = locStepIterator->second;
	}
	locParticleComboBlueprint->Prepend_ParticleComboBlueprintStep(locParticleComboBlueprintStep);

	--locStepIndex;
	if(dDebugLevel > 10)
		cout << "handle end: new step index, #steps = " << locStepIndex << ", " << locReaction->Get_NumReactionSteps() << endl;
	if(locStepIndex != -1)
	{
		// did not complete the chain yet
		locParticleIndex = 0;

		locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
		locParticleComboBlueprintStep->Set_ReactionStep(locReaction->Get_ReactionStep(locStepIndex));
		locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(locInitialParticleStepFromIndexMap[locStepIndex]);
		return true;
	}

	if(dDebugLevel > 10)
		cout << "save combo" << endl;
	_data.push_back(locParticleComboBlueprint);

	//register steps so they won't accidentally be recycled later, and so that they can be 
	for(size_t loc_i = 0; loc_i < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_i)
	{
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i);
		if(dSavedBlueprintSteps.find(locParticleComboBlueprintStep) != dSavedBlueprintSteps.end())
			continue;
		dSavedBlueprintSteps.insert(locParticleComboBlueprintStep);
		dBlueprintStepMap[*locParticleComboBlueprintStep] = const_cast<DParticleComboBlueprintStep*>(locParticleComboBlueprintStep);
	}

	locParticleComboBlueprint = new DParticleComboBlueprint(*locParticleComboBlueprint); //clone so don't alter saved object
	locParticleComboBlueprintStep = NULL;
	if(!Handle_Decursion(locParticleComboBlueprint, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleIndex, locStepIndex, locParticleComboBlueprintStep))
		return false;
	return true;
}

bool DParticleComboBlueprint_factory::Handle_Decursion(DParticleComboBlueprint* locParticleComboBlueprint, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque, int& locParticleIndex, int& locStepIndex, DParticleComboBlueprintStep*& locParticleComboBlueprintStep)
{
	do
	{
		if(dDebugLevel > 50)
			cout << "decursion: step, particle indices = " << locStepIndex << ", " << locParticleIndex << endl;
		if(locStepIndex == -1) //just saved a blueprint
		{
			if(dDebugLevel > 50)
				cout << "just saved blueprint" << endl;
			++locStepIndex;

			//this looks like it loses the pointer during the pop step, but it doesn't: it was just saved in dSavedBlueprintSteps
			locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
			*locParticleComboBlueprintStep = *(locParticleComboBlueprint->Pop_ParticleComboBlueprintStep());
			locParticleIndex = locResumeAtIndexDeque[locStepIndex].size() - 1;
			dCurrentComboSourceObjects.erase(locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject());
			if(dDebugLevel > 50)
				cout << "step index, particle index, resume at, #possible = " << locStepIndex << ", " << locParticleIndex << ", " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
			continue;
		}
		else if(locParticleIndex == int(locNumPossibilitiesDeque[locStepIndex].size())) //end of a step: step was found to be a duplicate OR failed a mass cut
		{
			if(dDebugLevel > 50)
				cout << "failed at end of a step" << endl;
			--locParticleIndex;
			dCurrentComboSourceObjects.erase(locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject());
			if(dDebugLevel > 50)
				cout << "step index, particle index, resume at, #possible = " << locStepIndex << ", " << locParticleIndex << ", " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
			continue;
		}

		//locParticleIndex will represent the particle it failed to find a combo for
		locResumeAtIndexDeque[locStepIndex][locParticleIndex] = 0;
		--locParticleIndex;
		if(locParticleIndex >= 0) //else is initial particle (will pop the entire step)
			dCurrentComboSourceObjects.erase(locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject());
		else
		{
			if(dDebugLevel > 50)
				cout << "pop the step" << endl;
			++locStepIndex;
			if(locStepIndex >= int(locResumeAtIndexDeque.size()))
				return false; //end of finding all blueprints

			const DParticleComboBlueprintStep* locPoppedStep = locParticleComboBlueprint->Pop_ParticleComboBlueprintStep();
			if(dSavedBlueprintSteps.find(locPoppedStep) == dSavedBlueprintSteps.end())
			{
				//popped step is not in a saved object: just use it, and recycle the step object we were just editing
				dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep);
				locParticleComboBlueprintStep = const_cast<DParticleComboBlueprintStep*>(locPoppedStep);
			}
			else //popped step is in a saved object: cannot edit it: use the current object
				*locParticleComboBlueprintStep = *locPoppedStep;

			locParticleIndex = locResumeAtIndexDeque[locStepIndex].size() - 1;
			dCurrentComboSourceObjects.erase(locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject());
		}
		if(dDebugLevel > 50)
			cout << "resume at, #possible = " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
	}
	while(locResumeAtIndexDeque[locStepIndex][locParticleIndex] == locNumPossibilitiesDeque[locStepIndex][locParticleIndex]);
	return true;
}

bool DParticleComboBlueprint_factory::Check_IfDuplicateStepCombo(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, int locStepIndex, deque<deque<int> >& locResumeAtIndexDeque, const deque<deque<int> >& locNumPossibilitiesDeque) const
{
#ifdef VTRACE
	VT_TRACER("DParticleComboBlueprint_factory::Check_IfDuplicateStepCombo()");
#endif
	//note that final particle ids could be rearranged in a different order
	map<Particle_t, unsigned int> locParticleTypeCount_CurrentStep;
	bool locIsAParticleDetected = false;
	for(size_t loc_i = 0; loc_i < locCurrentStep->Get_NumFinalParticleSourceObjects(); ++loc_i)
	{
		Particle_t locPID = locCurrentStep->Get_FinalParticleID(loc_i);
		if(locParticleTypeCount_CurrentStep.find(locPID) == locParticleTypeCount_CurrentStep.end())
			locParticleTypeCount_CurrentStep[locPID] = 1;
		else
			++(locParticleTypeCount_CurrentStep[locPID]);
		if(locCurrentStep->Is_FinalParticleDetected(loc_i))
			locIsAParticleDetected = true;
	}
	if(!locIsAParticleDetected)
		return false; //dupes of this sort only occur when dealing with at least some detected particles


	//if all of the particle types in this step are identical to all of the particle types in a previously-done step (and none of either are missing), regardless of the order they are listed:
		//the locResumeAtIndex system of the new step's particles must be >= the locResumeAtIndex system of the particles in the LAST OF the matching previously-done steps (e.g. could be 3x pi0s)
			//in other words: if 2x pi0 -> g, g; CD, AB ok and AB, CD bad; but also BC, AD ok (BC system occurs after AD system)
			//this works regardless of how many particles or what their PIDs are, or if it's a long decay chain or not

	//search future steps for a match (identical particles, none missing, may be in a different order)
	int locNumFutureSteps = locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps();
	for(int loc_i = 0; loc_i < locNumFutureSteps; ++loc_i)
	{
		int locFutureStepIndex = locResumeAtIndexDeque.size() - locNumFutureSteps + loc_i;
		const DParticleComboBlueprintStep* locFutureStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i);
		if(!Check_IfStepsAreIdentical(locParticleComboBlueprint, locCurrentStep, locFutureStep))
			continue;

		//step is a match
		if(dDebugLevel > 10)
			cout << "step at index = " << locStepIndex << " matches a future step at index " << locFutureStepIndex << endl;

		//the resume-at index of the first detected particle listed in the future step MUST be less than that of the first particle listed in the current step that has the same pid
			//else the current particle combo was already previously used for the future step in a previous combo
			//there must be at least one detected particle because this was required earlier
		size_t locFirstDetectedParticleIndex = 0;
		for(size_t loc_j = 0; loc_j < locFutureStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			if(!locFutureStep->Is_FinalParticleDetected(loc_j))
				continue;
			locFirstDetectedParticleIndex = loc_j;
			break;
		}

		Particle_t locFutureFirstPID = locFutureStep->Get_FinalParticleID(locFirstDetectedParticleIndex);
		int locFutureResumeAtIndex = locResumeAtIndexDeque[locFutureStepIndex][locFirstDetectedParticleIndex];
		for(size_t loc_j = 0; loc_j < locCurrentStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			if(locCurrentStep->Get_FinalParticleID(loc_j) != locFutureFirstPID)
				continue;
			int locCurrentResumeAtIndex = locResumeAtIndexDeque[locStepIndex][loc_j];
			if(dDebugLevel > 10)
				cout << "future first pid, future resume index, current resume index = " << locFutureFirstPID << ", " << locFutureResumeAtIndex << ", " << locCurrentResumeAtIndex << endl;
			if(locCurrentResumeAtIndex < locFutureResumeAtIndex)
			{
				if(loc_j == (locCurrentStep->Get_NumFinalParticleSourceObjects() - 1))
				{
					//on the last particle index: advance resume-at index to the smallest-possible, non-duplicate value to speed up the search
					locResumeAtIndexDeque[locStepIndex][loc_j] = locFutureResumeAtIndex;
					if(dDebugLevel > 10)
						cout << "resume-at index updated to " << locFutureResumeAtIndex << endl;
				}
				if(dDebugLevel > 10)
					cout << "duplicate step, force abort" << endl;
				return true;
			}

			if(dDebugLevel > 10)
				cout << "combo is ok (step is duplicate, but combo is not)" << endl;
			return false; //else combo is ok (step is duplicate, but combo is not)
		}
	}

	if(dDebugLevel > 10)
		cout << "step at index = " << locStepIndex << " has no similar steps." << endl;

	return false;
}

bool DParticleComboBlueprint_factory::Check_IfStepsAreIdentical(const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locCurrentStep, const DParticleComboBlueprintStep* locPreviousStep) const
{
	if(locCurrentStep->Get_MissingParticleIndex() != -1)
		return false; //missing particle somewhere in the step: cannot be a duplicate (already have checked not more than one missing particle)
	if(locPreviousStep->Get_MissingParticleIndex() != -1)
		return false; //missing particle somewhere in the step: cannot be a duplicate (already have checked not more than one missing particle)

	if(locPreviousStep->Get_InitialParticleID() != locCurrentStep->Get_InitialParticleID())
		return false;
	if(locPreviousStep->Get_TargetParticleID() != locCurrentStep->Get_TargetParticleID())
		return false;

	if(locPreviousStep->Get_NumFinalParticleSourceObjects() != locCurrentStep->Get_NumFinalParticleSourceObjects())
		return false;

	//a step is identical only if the decay chains of it's final state particles are also identical
	//note that final particle ids could be rearranged in a different order
	map<Particle_t, unsigned int> locParticleTypeCount_CurrentStep;
	map<Particle_t, deque<int> > locDecayingParticleStepIndices_CurrentStep;
	for(size_t loc_i = 0; loc_i < locCurrentStep->Get_NumFinalParticleSourceObjects(); ++loc_i)
	{
		Particle_t locPID = locCurrentStep->Get_FinalParticleID(loc_i);
		if(locParticleTypeCount_CurrentStep.find(locPID) == locParticleTypeCount_CurrentStep.end())
			locParticleTypeCount_CurrentStep[locPID] = 1;
		else
			++(locParticleTypeCount_CurrentStep[locPID]);
		if(locCurrentStep->Is_FinalParticleDecaying(loc_i))
			locDecayingParticleStepIndices_CurrentStep[locPID].push_back(locCurrentStep->Get_DecayStepIndex(loc_i));
	}

	//compare against the previous step
	map<Particle_t, deque<int> > locDecayingParticleStepIndices_PreviousStep;
	for(size_t loc_j = 0; loc_j < locPreviousStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
	{
		Particle_t locPID = locPreviousStep->Get_FinalParticleID(loc_j);
		if(locParticleTypeCount_CurrentStep.find(locPID) == locParticleTypeCount_CurrentStep.end())
			return false;
		if(locParticleTypeCount_CurrentStep[locPID] == 1)
			locParticleTypeCount_CurrentStep.erase(locPID); //if another one, will fail .find() check
		else
			--(locParticleTypeCount_CurrentStep[locPID]);
		if(locPreviousStep->Is_FinalParticleDecaying(loc_j))
			locDecayingParticleStepIndices_PreviousStep[locPID].push_back(locPreviousStep->Get_DecayStepIndex(loc_j));
	}
	if(!locParticleTypeCount_CurrentStep.empty())
		return false;

	//all of the particle types are the same, now check the decays
	map<Particle_t, deque<int> >::iterator locIterator;
	for(locIterator = locDecayingParticleStepIndices_CurrentStep.begin(); locIterator != locDecayingParticleStepIndices_CurrentStep.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		deque<int>& locDecayStepIndices_Current = locIterator->second;
		deque<int> locDecayStepIndices_Previous = locDecayingParticleStepIndices_PreviousStep[locPID];
		deque<int>::iterator locDequeIterator;
		for(size_t loc_i = 0; loc_i < locDecayStepIndices_Current.size(); ++loc_i)
		{
			bool locMatchFoundFlag = false;
			for(locDequeIterator = locDecayStepIndices_Previous.begin(); locDequeIterator != locDecayStepIndices_Previous.end(); ++locDequeIterator)
			{
				const DParticleComboBlueprintStep* locNewCurrentStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(locDecayStepIndices_Current[loc_i]);
				const DParticleComboBlueprintStep* locNewPreviousStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(*locDequeIterator);
				if(!Check_IfStepsAreIdentical(locParticleComboBlueprint, locNewCurrentStep, locNewPreviousStep))
					continue;
				locMatchFoundFlag = true;
				locDecayStepIndices_Previous.erase(locDequeIterator);
				break;
			}
			if(!locMatchFoundFlag)
				return false; //no matches found
		}
	}

	//step is a match
	return true;
}

int DParticleComboBlueprint_factory::Grab_DecayingParticle(Particle_t locAnalysisPID, int& locResumeAtIndex, const DReaction* locReaction, int locStepIndex, int locParticleIndex)
{
	if(locResumeAtIndex >= 1)
	{
		if(dDebugLevel > 10)
			cout << ParticleType(locAnalysisPID) << " does not decay later in the reaction." << endl;
		return -2;
	}

	if(dDebugLevel > 10)
		cout << "check if " << ParticleType(locAnalysisPID) << " decays later in the reaction." << endl;
	if((locAnalysisPID == Gamma) || (locAnalysisPID == Electron) || (locAnalysisPID == Positron) || (locAnalysisPID == Proton) || (locAnalysisPID == AntiProton))
	{
		if(dDebugLevel > 10)
			cout << ParticleType(locAnalysisPID) << " does not decay later in the reaction." << endl;
		return -2; //these particles don't decay: don't search!
	}

	//check to see how many final state particles with this pid type there are before now
	size_t locPreviousPIDCount = 0;
	for(int loc_i = 0; loc_i <= locStepIndex; ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			if((loc_i == locStepIndex) && (int(loc_j) == locParticleIndex))
				break; //at the current particle: of the search
			if(locReactionStep->Get_FinalParticleID(loc_j) == locAnalysisPID)
				++locPreviousPIDCount;
		}
	}

	//now, find the (locPreviousPIDCount + 1)'th time where this pid is a decay parent
	size_t locStepPIDCount = 0;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		if(locReactionStep->Get_InitialParticleID() != locAnalysisPID)
			continue;
		++locStepPIDCount;
		if(locStepPIDCount <= locPreviousPIDCount)
			continue;
		locResumeAtIndex = 1;
		if(dDebugLevel > 10)
			cout << ParticleType(locAnalysisPID) << " decays later in the reaction, at step index " << loc_i << endl;
		return loc_i;
	}

	if(dDebugLevel > 10)
		cout << ParticleType(locAnalysisPID) << " does not decay later in the reaction." << endl;
	return -2;
}

const JObject* DParticleComboBlueprint_factory::Grab_DetectedTrack(const DReaction* locReaction, Particle_t locAnalysisPID, int& locResumeAtIndex, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative)
{
	int locAnalysisCharge = ParticleCharge(locAnalysisPID);
	if(dDebugLevel > 10)
		cout << "Grab_DetectedTrack: PID, Charge = " << ParticleType(locAnalysisPID) << ", " << locAnalysisCharge << endl;
	if(locAnalysisCharge == 0)
		return Choose_SourceObject(locReaction, locAnalysisPID, locNeutralShowerDeque, locResumeAtIndex);
	else if(locAnalysisCharge > 0)
		return Choose_SourceObject(locReaction, locAnalysisPID, locChargedTrackDeque_Positive, locResumeAtIndex);
	else
		return Choose_SourceObject(locReaction, locAnalysisPID, locChargedTrackDeque_Negative, locResumeAtIndex);
}

const JObject* DParticleComboBlueprint_factory::Choose_SourceObject(const DReaction* locReaction, Particle_t locAnalysisPID, deque<const JObject*>& locSourceObjects, int& locResumeAtIndex) const
{
	if(dDebugLevel > 10)
		cout << "Choose_SourceObject: resume at, #possible = " << locResumeAtIndex << ", " << locSourceObjects.size() << endl;
	if(locResumeAtIndex >= int(locSourceObjects.size()))
		return NULL;
	const JObject* locObject = NULL;

	do
	{
		locObject = locSourceObjects[locResumeAtIndex];
		++locResumeAtIndex;

		//make sure not used currently
		if(dCurrentComboSourceObjects.find(locObject) != dCurrentComboSourceObjects.end())
		{
			if(dDebugLevel > 20)
				cout << "Source object already in use for locResumeAtIndex = " << locResumeAtIndex - 1 << endl;
			continue;
		}

		const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locObject); //NULL if not charged

		//if charged, check to make sure the tracking FOM is OK (cut garbage tracks and wildly bad combos)
		if(locChargedTrack != NULL)
		{
			bool locWillReSwimFlag = false;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_ChargedHypothesisToUse(locChargedTrack, locAnalysisPID, locWillReSwimFlag);
			if(!Cut_TrackingFOM(locReaction, locChargedTrackHypothesis))
			{
				if(dDebugLevel > 20)
					cout << "Bad Tracking FOM" << endl;
				continue;
			}
			else if(locWillReSwimFlag)
			{
				if(!Cut_HasDetectorMatch(locReaction, locChargedTrackHypothesis))
				{
					if(dDebugLevel > 20)
						cout << "No Detector Match" << endl;
					continue;
				}
			}
		}

		//check to make sure the track momentum isn't too low (e.g. testing a 100 MeV pion to be a proton)
		bool locTrackMomentumTooLowFlag = false;
		pair<bool, double> locMinProtonMomentum = dMinProtonMomentum.first ? dMinProtonMomentum : locReaction->Get_MinProtonMomentum();
		if((locChargedTrack != NULL) && locMinProtonMomentum.first && (ParticleMass(locAnalysisPID) >= (ParticleMass(Proton) - 0.001)))
		{
			if(locChargedTrack->Get_Hypothesis(Proton) == NULL)
			{
				deque<pair<Particle_t, bool> > locPIDsToTry = dTrackTimeBasedFactory_Combo->Get_ParticleIDsToTry(locAnalysisPID);
				bool locFoundFlag = false;
				for(size_t loc_i = 0; loc_i < locPIDsToTry.size(); ++loc_i)
				{
					const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPIDsToTry[loc_i].first);
					if(locChargedTrackHypothesis == NULL)
						continue;
					locFoundFlag = true;

					if(dDebugLevel > 20)
						cout << "Proton candidate, momentum = " << locChargedTrackHypothesis->momentum().Mag() << endl;
					if(locChargedTrackHypothesis->momentum().Mag() < locMinProtonMomentum.second)
						locTrackMomentumTooLowFlag = true;
					break;
				}
				if(!locFoundFlag)
				{
					const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
					if(dDebugLevel > 20)
						cout << "Proton candidate, momentum = " << locChargedTrackHypothesis->momentum().Mag() << endl;
					if(locChargedTrackHypothesis->momentum().Mag() < locMinProtonMomentum.second)
						locTrackMomentumTooLowFlag = true;
				}
			}
		}
		if(locTrackMomentumTooLowFlag)
		{
			if(dDebugLevel > 20)
				cout << "Track momentum too low to be " << ParticleType(locAnalysisPID) << endl;
			continue; //probably reconstructed a low-momentum pion: can't possibly be a proton (would stop too soon)
		}

		return locObject;
	}
	while(locResumeAtIndex < int(locSourceObjects.size()));
	return NULL;
}

const DChargedTrackHypothesis* DParticleComboBlueprint_factory::Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locAnalysisPID, bool& locWillReSwimFlag) const
{
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locAnalysisPID);
	if(locChargedTrackHypothesis != NULL)
		return locChargedTrackHypothesis;

	//pid not found for this track: loop over other possible pids
	deque<pair<Particle_t, bool> > locPIDsToTry = dTrackTimeBasedFactory_Combo->Get_ParticleIDsToTry(locAnalysisPID);
	for(size_t loc_i = 0; loc_i < locPIDsToTry.size(); ++loc_i)
	{
		locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPIDsToTry[loc_i].first);
		if(locChargedTrackHypothesis == NULL)
			continue;
		locWillReSwimFlag = locPIDsToTry[loc_i].second;
		return locChargedTrackHypothesis;
	}

	//still none found, take the one with the best FOM
	locWillReSwimFlag = true;
	return locChargedTrack->Get_BestFOM();
}

bool DParticleComboBlueprint_factory::Cut_HasDetectorMatch(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locHasDetectorMatchFlag = dHasDetectorMatchFlag.first ? dHasDetectorMatchFlag : locReaction->Get_HasDetectorMatchFlag();
	if((!locHasDetectorMatchFlag.first) || (!locHasDetectorMatchFlag.second))
		return true;

	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
	return dDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased);
}

bool DParticleComboBlueprint_factory::Cut_TrackingFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locMinTrackingFOM = dMinTrackingFOM.first ? dMinTrackingFOM : locReaction->Get_MinTrackingFOM();
	if(!locMinTrackingFOM.first)
		return true;
	double locFOM = TMath::Prob(locChargedTrackHypothesis->dChiSq_Track, locChargedTrackHypothesis->dNDF_Track);
	return ((locChargedTrackHypothesis->dNDF_Track == 0) ? true : (locFOM >= locMinTrackingFOM.second));
}

DParticleComboBlueprintStep* DParticleComboBlueprint_factory::Get_ParticleComboBlueprintStepResource(void)
{
	DParticleComboBlueprintStep* locParticleComboBlueprintStep;
	if(dParticleComboBlueprintStepPool_Available.empty())
	{
		locParticleComboBlueprintStep = new DParticleComboBlueprintStep;
		dParticleComboBlueprintStepPool_All.push_back(locParticleComboBlueprintStep);
	}
	else
	{
		locParticleComboBlueprintStep = dParticleComboBlueprintStepPool_Available.back();
		locParticleComboBlueprintStep->Reset();
		dParticleComboBlueprintStepPool_Available.pop_back();
	}
	return locParticleComboBlueprintStep;
}

void DParticleComboBlueprint_factory::Reset_Pools(void)
{
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dParticleComboBlueprintStepPool_All.size() > MAX_DParticleComboBlueprintStepPoolSize){
		for(size_t loc_i = MAX_DParticleComboBlueprintStepPoolSize; loc_i < dParticleComboBlueprintStepPool_All.size(); ++loc_i)
			delete dParticleComboBlueprintStepPool_All[loc_i];
		dParticleComboBlueprintStepPool_All.resize(MAX_DParticleComboBlueprintStepPoolSize);
	}
	dParticleComboBlueprintStepPool_Available = dParticleComboBlueprintStepPool_All;
}

bool DParticleComboBlueprint_factory::Calc_FinalStateP4(size_t locTotalNumSteps, const DParticleComboBlueprint* locParticleComboBlueprint, const DParticleComboBlueprintStep* locNewParticleComboBlueprintStep, int locStepIndexToGrab, DLorentzVector& locFinalStateP4) const
{
	//The input locParticleComboBlueprint is under construction: it does not have all of the steps yet
	//locNewParticleComboBlueprintStep is the step that will be added next, IF it passes this invariant mass cut
	//This is a recursive function, so locStepIndexToGrab is the index for the step that should be used for the p4 calculation at this stage
	//However, because the locParticleComboBlueprint is under construction, locStepIndexToGrab must be converted into the index that the step is actually located at
	const DParticleComboBlueprintStep* locCurrentStep = locNewParticleComboBlueprintStep;
	if(locStepIndexToGrab != -1)
	{
		int locStepFromBackIndex = locTotalNumSteps - locStepIndexToGrab;
		int locActualStepIndex = locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps() - locStepFromBackIndex;
		locCurrentStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(locActualStepIndex);
		if(locCurrentStep == NULL)
			return false;
	}

	locFinalStateP4.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
	for(size_t loc_i = 0; loc_i < locCurrentStep->Get_NumFinalParticleSourceObjects(); ++loc_i)
	{
		if(locCurrentStep->Is_FinalParticleMissing(loc_i))
			return false;
		if(locCurrentStep->Is_FinalParticleDecaying(loc_i))
		{
			DLorentzVector locDecayP4;
			if(!Calc_FinalStateP4(locTotalNumSteps, locParticleComboBlueprint, locNewParticleComboBlueprintStep, locCurrentStep->Get_DecayStepIndex(loc_i), locDecayP4))
				return false;
			locFinalStateP4 += locDecayP4;
		}
		else if(locCurrentStep->Is_FinalParticleCharged(loc_i))
		{
			const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locCurrentStep->Get_FinalParticle_SourceObject(loc_i));
			Particle_t locPID = locCurrentStep->Get_FinalParticleID(loc_i);
			bool locWillReSwimFlag = false;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_ChargedHypothesisToUse(locChargedTrack, locPID, locWillReSwimFlag);
			DVector3 locMomentum(locChargedTrackHypothesis->momentum());
			locFinalStateP4 += DLorentzVector(locMomentum, sqrt(locMomentum.Mag2() + ParticleMass(locPID)*ParticleMass(locPID)));
		}
		else //neutral
		{
			const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locCurrentStep->Get_FinalParticle_SourceObject(loc_i));
			DVector3 locHitPoint = locNeutralShower->dSpacetimeVertex.Vect();
			DVector3 locMomentum(locHitPoint - dVertex->dSpacetimeVertex.Vect());
			Particle_t locPID = locCurrentStep->Get_FinalParticleID(loc_i);
			if(locPID != Gamma)
			{
				double locDeltaT = locNeutralShower->dSpacetimeVertex.T() - dVertex->dSpacetimeVertex.T();
				double locBeta = locMomentum.Mag()/(locDeltaT*29.9792458); //path length is locMomentum.Mag() (for now)
				if(locBeta >= 1.0)
					locBeta = 0.9999;
				if(locBeta < 0.0)
					locBeta = 0.0;
				double locGamma = 1.0/sqrt(1.0 - locBeta*locBeta);
				double locMass = ParticleMass(locPID);
				double locPMag = locGamma*locBeta*locMass;
				locMomentum.SetMag(locPMag);
				locFinalStateP4 += DLorentzVector(locMomentum, sqrt(locMass*locMass + locPMag*locPMag));
			}
			else
			{
				locMomentum.SetMag(locNeutralShower->dEnergy);
				locFinalStateP4 += DLorentzVector(locMomentum, locNeutralShower->dEnergy);
			}
		}
	}
	return true;
}

//------------------
// erun
//------------------
jerror_t DParticleComboBlueprint_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleComboBlueprint_factory::fini(void)
{
	for(size_t loc_i = 0; loc_i < dParticleComboBlueprintStepPool_All.size(); ++loc_i)
		delete dParticleComboBlueprintStepPool_All[loc_i];

	return NOERROR;
}

