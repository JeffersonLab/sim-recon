#include "DParticleComboBlueprint_factory.h"

//------------------
// init
//------------------
jerror_t DParticleComboBlueprint_factory::init(void)
{
	dDebugLevel = 0;
	MAX_DParticleComboBlueprintStepPoolSize = 40;
	dMinimumProtonMomentum = 0.3; //from MIN_PROTON_P defined in DTrackFitterKalmanSIMD (as of August 12, 2012)
	dVertexZCutFlag = true;
	dMinVertexZ = 45.0;
	dMaxVertexZ = 85.0;
	dMaximumNumTracks = -1;
	dMinChargedPIDFOM = 0.001; //set to < 0.0 to disable
	dMaxTrackingChiSqPerDF = -1.0; //set to < 0.0 to disable
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleComboBlueprint_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("COMBOBLUEPRINTS:DEBUGLEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("COMBOBLUEPRINTS:MINPROTONMOMENTUM", dMinimumProtonMomentum);
	gPARMS->SetDefaultParameter("COMBOBLUEPRINTS:MAXNUMTRACKS", dMaximumNumTracks);
	gPARMS->SetDefaultParameter("COMBO:VERTEXZCUTFLAG", dVertexZCutFlag);
	gPARMS->SetDefaultParameter("COMBO:MINVERTEXZ", dMinVertexZ);
	gPARMS->SetDefaultParameter("COMBO:MAXVERTEXZ", dMaxVertexZ);
	gPARMS->SetDefaultParameter("COMBO:MINCHARGEDPIDFOM", dMinChargedPIDFOM);
	gPARMS->SetDefaultParameter("COMBO:MAXTRACKINGCHISQPERDF", dMaxTrackingChiSqPerDF);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleComboBlueprint_factory::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	Reset_Pools();

	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<const DReaction*> locReactions;
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		
		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		locReactions.insert(locReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}

	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		Build_ParticleComboBlueprints(locEventLoop, locReactions[loc_i]);

	return NOERROR;
}

jerror_t DParticleComboBlueprint_factory::Build_ParticleComboBlueprints(JEventLoop* locEventLoop, const DReaction* locReaction)
{
	vector<const DChargedTrack*> locChargedTrackVector;
	locEventLoop->Get(locChargedTrackVector);
	vector<const DNeutralShower*> locNeutralShowerVector;
	locEventLoop->Get(locNeutralShowerVector);

	if(dMaximumNumTracks >= 0)
	{
		if(int(locChargedTrackVector.size()) > dMaximumNumTracks)
			return RESOURCE_UNAVAILABLE;
	}

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
	for(size_t loc_i = 0; loc_i < locChargedTrackVector.size(); ++loc_i)
	{
		locChargedTrack = locChargedTrackVector[loc_i];
		if(locChargedTrack->Get_Charge() > 0.0)
			locChargedTrackDeque_Positive.push_back(static_cast<const JObject*>(locChargedTrack));
		else
			locChargedTrackDeque_Negative.push_back(static_cast<const JObject*>(locChargedTrack));
	}

	if(dDebugLevel > 0)
		cout << "#+, #-, #0 particles = " << locChargedTrackDeque_Positive.size() << ", " << locChargedTrackDeque_Negative.size() << ", " << locNeutralShowerDeque.size() << endl;

	//set up combo loop
	deque<deque<int> > locResumeAtIndexDeque; //1st index is step, 2nd is particle (the initial particle, then the final particles)
	deque<deque<int> > locNumPossibilitiesDeque; //1st index is step, 2nd is particle (the initial particle, then the final particles)
	Setup_ComboLoop(locReaction, locNeutralShowerDeque.size(), locChargedTrackDeque_Positive.size(), locChargedTrackDeque_Negative.size(), locResumeAtIndexDeque, locNumPossibilitiesDeque);

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
	}

	//find the combos!!
	vector<DParticleComboBlueprint*> locParticleComboBlueprints;
	Find_Combos(locReaction, locNeutralShowerDeque, locChargedTrackDeque_Positive, locChargedTrackDeque_Negative, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleComboBlueprints);

	if(dDebugLevel > 10)
	{
		cout << "print pointers: " << endl;
		for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
		{
			cout << "COMBO " << loc_i << endl;
			for(size_t loc_j = 0; loc_j < locParticleComboBlueprints[loc_i]->Get_NumParticleComboBlueprintSteps(); ++loc_j)
			{
				cout << "Step " << loc_j << " pointers: ";
				for(size_t loc_k = 0; loc_k < locParticleComboBlueprints[loc_i]->Get_ParticleComboBlueprintStep(loc_j)->Get_NumFinalParticleSourceObjects(); ++loc_k)
					cout << locParticleComboBlueprints[loc_i]->Get_ParticleComboBlueprintStep(loc_j)->Get_FinalParticle_SourceObject(loc_k) << ", ";
				cout << endl;
			}
		}
	}

	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
		_data.push_back(locParticleComboBlueprints[loc_i]);

	return NOERROR;
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

void DParticleComboBlueprint_factory::Find_Combos(const DReaction* locReaction, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints)
{
	DParticleComboBlueprint* locParticleComboBlueprint = new DParticleComboBlueprint();
	locParticleComboBlueprint->Set_Reaction(locReaction);

	int locStepIndex = 0;
	int locParticleIndex = 0; //final = 0 -> (#final - 1)
	DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
	locParticleComboBlueprintStep->Set_ReactionStep(locReaction->Get_ReactionStep(0));
	locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(-1);

	do{
		if(dDebugLevel > 10)
			cout << "do loop: step & particle indices = " << locStepIndex << ", " << locParticleIndex << endl;
		if(locParticleIndex == int(locNumPossibilitiesDeque[locStepIndex].size()))
		{
			if(dDebugLevel > 10)
				cout << "end of reaction step" << endl;
			if(!Handle_EndOfReactionStep(locReaction, locParticleComboBlueprint, locParticleComboBlueprintStep, locStepIndex, locParticleIndex, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleComboBlueprints))
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

		//if two particles of same type in a step: locResumeAtIndex must always be >= the previous locResumeAtIndex (prevents duplicates) e.g. g, D -> p, p, pi-
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

		// check to see if this particle has a decay that is represented in a future step
		// e.g. on Lambda in g, p -> K+, Lambda; where a later step is Lambda -> p, pi-
		int locDecayStepIndex = Grab_DecayingParticle(locParticleComboBlueprint, locAnalysisPID, locResumeAtIndex, locReaction, locStepIndex, locParticleIndex, locParticleComboBlueprintStep);
		if(locDecayStepIndex >= 0)
		{
			if(dDebugLevel > 10)
				cout << "decaying particle" << endl;
			locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(NULL, locDecayStepIndex); //decaying
			++locParticleIndex;
			continue;
		}

		// else grab a detected track
		const JObject* locSourcObject = Grab_DetectedTrack(locParticleComboBlueprint, locAnalysisPID, locResumeAtIndex, locNeutralShowerDeque, locChargedTrackDeque_Positive, locChargedTrackDeque_Negative);
		if(locSourcObject == NULL)
		{
			if(dDebugLevel > 10)
				cout << "can't find detected particle" << endl;
			if(!Handle_Decursion(locParticleComboBlueprint, locResumeAtIndexDeque, locNumPossibilitiesDeque, locParticleIndex, locStepIndex, locParticleComboBlueprintStep))
				break;
			continue;
		}

		if(dDebugLevel > 10)
			cout << "detected track found, locResumeAtIndex now = " << locResumeAtIndex << endl;
		locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(locSourcObject, -2); //detected
		++locParticleIndex;
	}
	while(true);
	delete locParticleComboBlueprint; //delete the last, extra one
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

DParticleComboBlueprint* DParticleComboBlueprint_factory::Clone_ParticleComboBlueprint(const DParticleComboBlueprint* locParticleComboBlueprint)
{
	DParticleComboBlueprint* locNewParticleComboBlueprint = new DParticleComboBlueprint();
	*locNewParticleComboBlueprint = *locParticleComboBlueprint;
	return locNewParticleComboBlueprint;
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

bool DParticleComboBlueprint_factory::Handle_EndOfReactionStep(const DReaction* locReaction, DParticleComboBlueprint*& locParticleComboBlueprint, DParticleComboBlueprintStep*& locParticleComboBlueprintStep, int& locStepIndex, int& locParticleIndex, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque, vector<DParticleComboBlueprint*>& locParticleComboBlueprints)
{
	//end of step, advance to next step

	//first check to see if identical to a previous step; if so, just save the old step and recycle the current one
	bool locRecycleFlag = false;
	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprints[loc_i]->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			//could have two identical steps in the same reaction (e.g. pi0 decays)
			const DParticleComboBlueprintStep* locPreviousParticleComboBlueprintStep = locParticleComboBlueprints[loc_i]->Get_ParticleComboBlueprintStep(loc_j);
			if((*locPreviousParticleComboBlueprintStep) != (*locParticleComboBlueprintStep))
				continue; //else identical: just use the previously saved one and recycle the newly constructed one
			if(dDebugLevel > 10)
				cout << "step identical to previous one: copying & recycling." << endl;
			Recycle_ParticleComboBlueprintStep(locParticleComboBlueprintStep);
			locParticleComboBlueprintStep = NULL;
			locParticleComboBlueprint->Add_ParticleComboBlueprintStep(locPreviousParticleComboBlueprintStep);
			locRecycleFlag = true;
			break;
		}
		if(locRecycleFlag)
			break;
	}
	//if match not found, try searching previous dreactions too (may be identical reactions or similar enough)
	if(!locRecycleFlag)
	{
		for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < _data[loc_i]->Get_NumParticleComboBlueprintSteps(); ++loc_j)
			{
				const DParticleComboBlueprintStep* locPreviousParticleComboBlueprintStep = _data[loc_i]->Get_ParticleComboBlueprintStep(loc_j);
				if((*locPreviousParticleComboBlueprintStep) != (*locParticleComboBlueprintStep))
					continue; //else identical: just use the previously saved one and recycle the newly constructed one
				if(dDebugLevel > 10)
					cout << "step identical to one from previous DReaction: copying & recycling." << endl;
				Recycle_ParticleComboBlueprintStep(locParticleComboBlueprintStep);
				locParticleComboBlueprintStep = NULL;
				locParticleComboBlueprint->Add_ParticleComboBlueprintStep(locPreviousParticleComboBlueprintStep);
				locRecycleFlag = true;
				break;
			}
			if(locRecycleFlag)
				break;
		}
	}
	if(!locRecycleFlag)
		locParticleComboBlueprint->Add_ParticleComboBlueprintStep(locParticleComboBlueprintStep);

	++locStepIndex;
	if(dDebugLevel > 10)
		cout << "handle end: new step index, #steps = " << locStepIndex << ", " << locReaction->Get_NumReactionSteps() << endl;
	if(locStepIndex != int(locReaction->Get_NumReactionSteps()))
	{
		// did not complete the chain yet
		locParticleIndex = 0;

		locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
		locParticleComboBlueprintStep->Set_ReactionStep(locReaction->Get_ReactionStep(locStepIndex));
		//loop back through previous steps, see where the initial particle of the next (not-just-finished) step decayed from
		bool locMatchFoundFlag = false;
		for(int loc_i = 0; loc_i < locStepIndex; ++loc_i)
		{
			for(int loc_j = 0; loc_j < int(locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i)->Get_NumFinalParticleSourceObjects()); ++loc_j)
			{
				if(dDebugLevel > 10)
					cout << "previous step index, previous step final particle index, dDecayStepIndex, locStepIndex = " << loc_i << ", " << loc_j << ", " << locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i)->Get_DecayStepIndex(loc_j) << ", " << locStepIndex << endl;
				if(locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i)->Get_DecayStepIndex(loc_j) == locStepIndex)
				{
					locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(loc_i);
					if(dDebugLevel > 10)
						cout << "initial particle (" << ParticleType(locParticleComboBlueprintStep->Get_InitialParticleID()) << ") found in previous step index: " << loc_i << endl;
					locMatchFoundFlag = true;
					break;
				}
			}
			if(locMatchFoundFlag)
				break;
		}
		if(!locMatchFoundFlag)
			return false; //reaction setup incorrectly: just break

		return true;
	}

	//constructed, handle decursion
	locParticleComboBlueprints.push_back(locParticleComboBlueprint);
	locParticleComboBlueprint = Clone_ParticleComboBlueprint(locParticleComboBlueprint);
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
		if(locStepIndex == int(locResumeAtIndexDeque.size())) //just saved a test reaction
		{
			if(dDebugLevel > 50)
				cout << "saved test reaction" << endl;
			--locStepIndex;
			locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
			*locParticleComboBlueprintStep = *(locParticleComboBlueprint->Pop_ParticleComboBlueprintStep());
			locParticleIndex = locResumeAtIndexDeque[locStepIndex].size() - 1;
			locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject();
			if(dDebugLevel > 50)
				cout << "step index, particle index, resume at, #possible = " << locStepIndex << ", " << locParticleIndex << ", " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
			continue;
		}
		else if(locParticleIndex == int(locNumPossibilitiesDeque[locStepIndex].size())) //end of a step: step was found to be a duplicate
		{
			//is this even possible anymore?
			if(dDebugLevel > 50)
				cout << "end of a step" << endl;
			--locParticleIndex;
			locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject();
			if(dDebugLevel > 50)
				cout << "step index, particle index, resume at, #possible = " << locStepIndex << ", " << locParticleIndex << ", " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
			continue;
		}

		//locParticleIndex will represent the particle it failed to find a combo for
		locResumeAtIndexDeque[locStepIndex][locParticleIndex] = 0;
		--locParticleIndex;
		if(locParticleIndex >= 0) //else is initial particle (will pop the entire step)
			locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject();
		else
		{
			if(dDebugLevel > 50)
				cout << "pop the step" << endl;
			--locStepIndex;
			if(locStepIndex < 0)
				return false;
			dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep); //recycle step!
			locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
			*locParticleComboBlueprintStep = *(locParticleComboBlueprint->Pop_ParticleComboBlueprintStep());
			locParticleIndex = locResumeAtIndexDeque[locStepIndex].size() - 1;
			locParticleComboBlueprintStep->Pop_FinalParticle_SourceObject();
		}
		if(dDebugLevel > 50)
			cout << "resume at, #possible = " << locResumeAtIndexDeque[locStepIndex][locParticleIndex] << ", " << locNumPossibilitiesDeque[locStepIndex][locParticleIndex] << endl;
	}
	while(locResumeAtIndexDeque[locStepIndex][locParticleIndex] == locNumPossibilitiesDeque[locStepIndex][locParticleIndex]);
	return true;
}

void DParticleComboBlueprint_factory::Setup_ComboLoop(const DReaction* locReaction, int locNumDetectedNeutralParticles, int locNumDetectedPositiveParticles, int locNumDetectedNegativeParticles, deque<deque<int> >& locResumeAtIndexDeque, deque<deque<int> >& locNumPossibilitiesDeque)
{
	//setup locResumeAtIndexDeque, & locNumPossibilitiesDeque
	Particle_t locAnalysisPID;
	int locMissingParticleIndex;
	unsigned int locNumSteps = locReaction->Get_NumReactionSteps();
	locResumeAtIndexDeque.clear();
	locResumeAtIndexDeque.clear();
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

			//check to see if it's decay is represented in a future step
			locTempDeque[loc_j] = 0;
			for(size_t loc_k = loc_i + 1; loc_k < locNumSteps; ++loc_k)
			{
				if(locReaction->Get_ReactionStep(loc_k)->Get_InitialParticleID() == locAnalysisPID)
					++(locTempDeque[loc_j]);
			}
			if(locTempDeque[loc_j] > 0)
				continue; //does decay, got the #combos

			//else use detected particles
			if(ParticleCharge(locAnalysisPID) == 0)
				locTempDeque[loc_j] = locNumDetectedNeutralParticles;
			else if(ParticleCharge(locAnalysisPID) == 1)
				locTempDeque[loc_j] = locNumDetectedPositiveParticles;
			else
				locTempDeque[loc_j] = locNumDetectedNegativeParticles;
		}
		locNumPossibilitiesDeque.push_back(locTempDeque);
	}
}

int DParticleComboBlueprint_factory::Grab_DecayingParticle(DParticleComboBlueprint* locParticleComboBlueprint, Particle_t locAnalysisPID, int& locResumeAtIndex, const DReaction* locReaction, int locStepIndex, int locParticleIndex, DParticleComboBlueprintStep* locParticleComboBlueprintStep)
{
	if(dDebugLevel > 10)
		cout << "check if " << ParticleType(locAnalysisPID) << " decays later in the reaction." << endl;
	if((locAnalysisPID == Gamma) || (locAnalysisPID == Electron) || (locAnalysisPID == Positron) || (locAnalysisPID == Proton) || (locAnalysisPID == AntiProton))
	{
		if(dDebugLevel > 10)
			cout << ParticleType(locAnalysisPID) << " does not decay later in the reaction." << endl;
		return -2; //these particles don't decay: don't search!
	}

	int locCurrentIndex = -1;
	const DReactionStep* locReactionStep;
	for(size_t loc_i = locStepIndex + 1; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		locReactionStep = locReaction->Get_ReactionStep(loc_i);
		if(locReactionStep->Get_InitialParticleID() != locAnalysisPID)
			continue;

		//have a match, but skip it if not at least at the correct resume index
		++locCurrentIndex;
		if(dDebugLevel > 10)
			cout << "match: current, resume at indices = " << locCurrentIndex << ", " << locResumeAtIndex << endl;
		if(locCurrentIndex < locResumeAtIndex)
			continue;

		//check to make sure not already used previously
		bool locPreviouslyUsedFlag = false;
		for(int loc_j = 0; loc_j < locStepIndex; ++loc_j)
		{
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j)->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				if(dDebugLevel > 15)
					cout << "i, j, k, dDecayStepIndex = " << loc_i << ", " << loc_j << ", " << loc_k << ", " << locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j)->Get_DecayStepIndex(loc_k) << endl;
				if(locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j)->Get_DecayStepIndex(loc_k) == int(loc_i))
				{
					if(dDebugLevel > 15)
						cout << "used previously" << endl;
					locPreviouslyUsedFlag = true;
					break;
				}
			}
			if(locPreviouslyUsedFlag)
				break;
		}
		if(locPreviouslyUsedFlag)
			continue; //this one has been previously used: continue;

		//check to see if it's been used previously in the current step also:
		for(int loc_j = 0; loc_j < locParticleIndex; ++loc_j)
		{
			if(dDebugLevel > 15)
				cout << "i, j, dDecayStepIndex = " << loc_i << ", " << loc_j << ", " << locParticleComboBlueprintStep->Get_DecayStepIndex(loc_j) << endl;
			if(locParticleComboBlueprintStep->Get_DecayStepIndex(loc_j) == int(loc_i))
			{
				if(dDebugLevel > 15)
					cout << "used previously" << endl;
				locPreviouslyUsedFlag = true;
				break;
			}
		}
		if(locPreviouslyUsedFlag)
			continue; //this one has been previously used: continue;

		//else it's good!!
		locResumeAtIndex = locCurrentIndex + 1;
		if(dDebugLevel > 10)
			cout << ParticleType(locAnalysisPID) << " decays later in the reaction, at step index " << loc_i << endl;

		return loc_i;
	}

	if(dDebugLevel > 10)
		cout << ParticleType(locAnalysisPID) << " does not decay later in the reaction." << endl;
	return -2;
}

const JObject* DParticleComboBlueprint_factory::Grab_DetectedTrack(DParticleComboBlueprint* locParticleComboBlueprint, Particle_t locAnalysisPID, int& locResumeAtIndex, deque<const JObject*>& locNeutralShowerDeque, deque<const JObject*>& locChargedTrackDeque_Positive, deque<const JObject*>& locChargedTrackDeque_Negative)
{
	int locAnalysisCharge = ParticleCharge(locAnalysisPID);
	if(dDebugLevel > 10)
		cout << "Grab_DetectedTrack: PID, Charge = " << ParticleType(locAnalysisPID) << ", " << locAnalysisCharge << endl;
	if(locAnalysisCharge == 0)
		return Choose_SourceObject(locAnalysisPID, locParticleComboBlueprint, locNeutralShowerDeque, locResumeAtIndex);
	else if(locAnalysisCharge > 0)
		return Choose_SourceObject(locAnalysisPID, locParticleComboBlueprint, locChargedTrackDeque_Positive, locResumeAtIndex);
	else
		return Choose_SourceObject(locAnalysisPID, locParticleComboBlueprint, locChargedTrackDeque_Negative, locResumeAtIndex);
}

const JObject* DParticleComboBlueprint_factory::Choose_SourceObject(Particle_t locAnalysisPID, DParticleComboBlueprint* locParticleComboBlueprint, deque<const JObject*>& locSourceObjects, int& locResumeAtIndex) const
{
	if(dDebugLevel > 10)
		cout << "Choose_SourceObject: resume at, #possible = " << locResumeAtIndex << ", " << locSourceObjects.size() << endl;
	if(locResumeAtIndex >= int(locSourceObjects.size()))
		return NULL;
	const JObject* locObject;

	const DParticleComboBlueprintStep* locParticleComboBlueprintStep;
	do
	{
		locObject = locSourceObjects[locResumeAtIndex];
		//make sure not used currently
		bool locTrackInUseFlag = false;
		for(size_t loc_i = 0; loc_i < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_i)
		{
			locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_i);
			for(size_t loc_j = 0; loc_j < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
			{
				if(locObject == locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_j))
				{
					if(dDebugLevel > 20)
						cout << "Source object already in use for locResumeAtIndex = " << locResumeAtIndex << endl;
					locTrackInUseFlag = true;
					break;
				}
			}
			if(locTrackInUseFlag)
				break;
		}
		++locResumeAtIndex;
		if(locTrackInUseFlag)
			continue;

		const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locObject); //NULL if not charged


		//if charged, check to make sure the vertex-z, PID FOM, and tracking chisq are OK (cut garbage tracks and wildly bad combos)
		if(locChargedTrack != NULL)
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locAnalysisPID);
			if(locChargedTrackHypothesis != NULL)
			{
				if((!Cut_VertexZ(locChargedTrackHypothesis)) || (!Cut_PIDFOM(locChargedTrackHypothesis)) || (!Cut_TrackingChiSqPerDF(locChargedTrackHypothesis)))
				{
					if(dDebugLevel > 20)
						cout << "Bad vertex-z, PID FOM, or Tracking chisq/df: " << locChargedTrackHypothesis->position().Z() << ", " << locChargedTrackHypothesis->dFOM << endl;
					continue;
				}
			}
			//pid not present: will need to swim: make sure that at least one of the hypotheses has a good vertex-z
				//the specific hypothesis used is determined by DTrackTimeBased_factory_Reaction, don't guess it here!!
					//will have to also cut later: in DParticleCombo_factory_PreKinFit
			bool locAllBadVertexZFlag = true;
			for(size_t loc_i = 0; loc_i < locChargedTrack->dChargedTrackHypotheses.size(); ++loc_i)
			{
				if(Cut_VertexZ(locChargedTrack->dChargedTrackHypotheses[loc_i]))
				{
					locAllBadVertexZFlag = false;
					break;
				}
			}
			if(locAllBadVertexZFlag)
			{
				if(dDebugLevel > 20)
					cout << "Bad vertex-z: " << locChargedTrack->Get_BestFOM()->position().Z() << endl;
				continue;
			}
		}

		//check to make sure the track momentum isn't too low (e.g. testing a 100 MeV pion to be a proton)
		bool locTrackMomentumTooLowFlag = false;
		if((locChargedTrack != NULL) && (ParticleMass(locAnalysisPID) >= ParticleMass(Proton)))
		{
			if(locChargedTrack->Get_Hypothesis(Proton) == NULL)
			{
				if(locChargedTrack->Get_Hypothesis(KPlus) != NULL)
				{
					if(dDebugLevel > 20)
						cout << "Proton candidate: K+ momentum = " << locChargedTrack->Get_Hypothesis(KPlus)->momentum().Mag() << endl;
					if(locChargedTrack->Get_Hypothesis(KPlus)->momentum().Mag() < dMinimumProtonMomentum)
						locTrackMomentumTooLowFlag = true;
				}
				else
				{
					if(dDebugLevel > 20)
						cout << "Proton candidate: Best FOM momentum = " << locChargedTrack->Get_BestFOM()->momentum().Mag() << endl;
					if(locChargedTrack->Get_BestFOM()->momentum().Mag() < dMinimumProtonMomentum)
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


bool DParticleComboBlueprint_factory::Cut_TrackingChiSqPerDF(const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	if(dMaxTrackingChiSqPerDF < 0.0)
		return true;
	if(locChargedTrackHypothesis->dNDF_Track == 0)
		return true;
	double locFOM = locChargedTrackHypothesis->dChiSq_Track/((double)(locChargedTrackHypothesis->dNDF_Track));
	return (locFOM <= dMaxTrackingChiSqPerDF);
}

bool DParticleComboBlueprint_factory::Cut_PIDFOM(const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	if(dMinChargedPIDFOM < 0.0)
		return true;
	return ((locChargedTrackHypothesis->dNDF == 0) ? true : (locChargedTrackHypothesis->dFOM >= dMinChargedPIDFOM));
}

bool DParticleComboBlueprint_factory::Cut_VertexZ(const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	if(!dVertexZCutFlag)
		return true;
	double locVertexZ = locChargedTrackHypothesis->position().Z();
	return ((locVertexZ >= dMinVertexZ) && (locVertexZ <= dMaxVertexZ));
}

