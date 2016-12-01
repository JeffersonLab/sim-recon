#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DParticleCombo_factory.h"

//------------------
// init
//------------------
jerror_t DParticleCombo_factory::init(void)
{
	MAX_DParticleComboStepPoolSize = 200;
	MAX_DKinematicDataPoolSize = 40;

	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
	//All combos that fail the kinematic fit, or aren't fit at all, are identical to the pre-combo versions.  
		//For these combos, just save the pointers to the previous objects in _data.  
	SetFactoryFlag(NOT_OBJECT_OWNER);
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleCombo_factory::brun(jana::JEventLoop* locEventLoop, int32_t runnumber)
{
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);

	// Get # of DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	size_t locNumReactions = 0;
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
		locNumReactions += locReactionsSubset.size();
	}

	MAX_DParticleComboStepPoolSize = 100*locNumReactions;
	MAX_DKinematicDataPoolSize = 100*locNumReactions;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleCombo_factory::evnt(JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DParticleCombo_factory::evnt()");
#endif

	Reset_Data();
	Reset_Pools();

 	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector, "PreKinFit");

	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

	set<const DParticleCombo*> locParticleCombos_FailedKinFit; //fill with all eligible, will erase as they are validated
	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		deque<const DParticleCombo*> locSurvivedParticleCombos;
		locAnalysisResultsVector[loc_i]->Get_PassedParticleCombos(locSurvivedParticleCombos);

		const DReaction* locReaction = locAnalysisResultsVector[loc_i]->Get_Reaction();
		if(locReaction->Get_KinFitType() == d_NoFit) //kinfit not requested, just clone the original objects
		{
			for(size_t loc_j = 0; loc_j < locSurvivedParticleCombos.size(); ++loc_j)
				_data.push_back(const_cast<DParticleCombo*>(locSurvivedParticleCombos[loc_j]));
			continue;
		}
		for(size_t loc_j = 0; loc_j < locSurvivedParticleCombos.size(); ++loc_j)
			locParticleCombos_FailedKinFit.insert(locSurvivedParticleCombos[loc_j]);
	}

	if(locKinFitResultsVector.empty())
	{
		//all failed! (or none requested) save 'em and bail
		set<const DParticleCombo*>::iterator locIterator = locParticleCombos_FailedKinFit.begin();
		for(; locIterator != locParticleCombos_FailedKinFit.end(); ++locIterator)
			_data.push_back(const_cast<DParticleCombo*>(*locIterator));
		return NOERROR; //kinfit not requested (or all kinfits failed): done
	}

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses, "KinFit");

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses, "KinFit");

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "KinFit");

	//map from old combo to new combo
	map<const DParticleCombo*, DParticleCombo*> locNewParticleComboMap;

	//loop over kinfit results
	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		map<const DParticleCombo*, const DKinFitChain*> locParticleComboMap;
		locKinFitResultsVector[loc_i]->Get_ParticleComboMap(locParticleComboMap); //each result may be for several combos (duplicate results)

		map<const DParticleCombo*, const DKinFitChain*>::iterator locComboIterator = locParticleComboMap.begin();
		for(; locComboIterator != locParticleComboMap.end(); ++locComboIterator)
		{
			const DParticleCombo* locParticleCombo = locComboIterator->first;
			const DKinFitChain* locKinFitChain = locComboIterator->second;

			locParticleCombos_FailedKinFit.erase(locParticleCombo); //kinfit successful, don't copy pointer later

			DParticleCombo* locNewParticleCombo = new DParticleCombo();
			dCreatedParticleCombos.push_back(locNewParticleCombo);
			locNewParticleCombo->Set_Reaction(locParticleCombo->Get_Reaction());
			locNewParticleCombo->Set_KinFitResults(locKinFitResultsVector[loc_i]);
			locNewParticleCombo->Set_EventRFBunch(locParticleCombo->Get_EventRFBunch());

			vector<const DParticleComboBlueprint*> locParticleComboBlueprints_Associated;
			locParticleCombo->GetT(locParticleComboBlueprints_Associated);
			locNewParticleCombo->AddAssociatedObject(locParticleComboBlueprints_Associated[0]);
			locNewParticleCombo->AddAssociatedObject(locParticleCombo);

			//search to see if the results were the same as a previous combo. if so, copy them:
				//note that it is possible for two combos with different steps to have the same kinfit results (e.g. one has an omega or phi, and the other doesn't)
			const DParticleCombo* locPreviousNewParticleCombo = Check_IsDuplicateCombo(locNewParticleComboMap, locParticleCombo);
			if(locPreviousNewParticleCombo != NULL)
			{
				//everything is identical: copy results
				for(size_t loc_k = 0; loc_k < locPreviousNewParticleCombo->Get_NumParticleComboSteps(); ++loc_k)
					locNewParticleCombo->Add_ParticleComboStep(locPreviousNewParticleCombo->Get_ParticleComboStep(loc_k));

				locNewParticleComboMap[locParticleCombo] = locNewParticleCombo;
				_data.push_back(locNewParticleCombo);
				continue;
			}

			//unique combo: build steps
			for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
			{
				const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
				DParticleComboStep* locNewParticleComboStep = Get_ParticleComboStepResource();

				locNewParticleComboStep->Set_MeasuredParticleComboStep(locParticleComboStep);
				locNewParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboStep->Get_ParticleComboBlueprintStep());
				locNewParticleComboStep->Set_SpacetimeVertex(locParticleComboStep->Get_SpacetimeVertex()); //overridden if kinematic fit

				//INITIAL PARTICLE & SPACETIME VERTEX
				if(locParticleComboStep->Is_InitialParticleDetected()) //set beam photon
				{
					bool locParticleFoundFlag = false;
					for(size_t loc_k = 0; loc_k < locBeamPhotons.size(); ++loc_k)
					{
						if(!locBeamPhotons[loc_k]->IsAssociated(locParticleCombo))
							continue; // beam photon created for a different particle combo
						locNewParticleComboStep->Set_InitialParticle(locBeamPhotons[loc_k]);
						locParticleFoundFlag = true;
						break;
					}
					if(!locParticleFoundFlag) //beam not used in kinfit, re-set original
						locNewParticleComboStep->Set_InitialParticle(locParticleComboStep->Get_InitialParticle());
				}
				else //decaying particle! //set here for initial state, and in previous step for final state
					Set_DecayingParticles(locNewParticleCombo, locParticleCombo, loc_j, locNewParticleComboStep, locKinFitChain, locKinFitResultsVector[loc_i]);

				//TARGET PARTICLE
				if(locParticleComboStep->Is_TargetPresent())
				{
					Particle_t locTargetPID = locParticleComboStep->Get_TargetParticleID();
					set<DKinFitParticle*> locTargetParticles = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticles(d_TargetParticle);
					if(!locTargetParticles.empty())
					{
						DKinematicData* locNewKinematicData = Build_KinematicData(locTargetPID, *locTargetParticles.begin());
						locNewParticleComboStep->Set_TargetParticle(locNewKinematicData);
					}
					else //not used in kinfit, re-set the original
						locNewParticleComboStep->Set_TargetParticle(locParticleComboStep->Get_TargetParticle());
				}
				else
					locNewParticleComboStep->Set_TargetParticle(NULL);

				//FINAL PARTICLES
				for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
				{
					const DKinematicData* locKinematicData_Measured = locParticleComboStep->Get_FinalParticle_Measured(loc_k);
					Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_k);
					if(locParticleComboStep->Is_FinalParticleMissing(loc_k)) //missing!
					{
						set<DKinFitParticle*> locMissingParticles = locKinFitResultsVector[loc_i]->Get_OutputKinFitParticles(d_MissingParticle);
						if(!locMissingParticles.empty())
						{
							DKinematicData* locNewKinematicData = Build_KinematicData(locPID, *locMissingParticles.begin());
							locNewParticleComboStep->Add_FinalParticle(locNewKinematicData);
						}
						else //not used in kinfit: do not create: NULL
							locNewParticleComboStep->Add_FinalParticle(NULL);
					}
					else if(locParticleComboStep->Is_FinalParticleDecaying(loc_k)) //decaying
						locNewParticleComboStep->Add_FinalParticle(NULL); //is set later, when it's in the initial state
					else if(locParticleComboStep->Is_FinalParticleNeutral(loc_k)) //neutral
						locNewParticleComboStep->Add_FinalParticle(Get_NeutralHypothesis(locParticleCombo, locNeutralParticleHypotheses, locKinematicData_Measured));
					else //charged
						locNewParticleComboStep->Add_FinalParticle(Get_ChargedHypothesis(locParticleCombo, locChargedTrackHypotheses, locKinematicData_Measured));
				}

				Set_SpacetimeVertex(locNewParticleCombo, locNewParticleComboStep, loc_j, locKinFitResultsVector[loc_i], locKinFitChain);
				locNewParticleCombo->Add_ParticleComboStep(locNewParticleComboStep);
			}

			locNewParticleComboMap[locParticleCombo] = locNewParticleCombo;
			_data.push_back(locNewParticleCombo);
		}
	}

	//directly save all combos for which the kinfits failed
	set<const DParticleCombo*>::iterator locIterator = locParticleCombos_FailedKinFit.begin();
	for(; locIterator != locParticleCombos_FailedKinFit.end(); ++locIterator)
		_data.push_back(const_cast<DParticleCombo*>(*locIterator));

	return NOERROR;
}

const DParticleCombo* DParticleCombo_factory::Check_IsDuplicateCombo(const map<const DParticleCombo*, DParticleCombo*>& locNewParticleComboMap, const DParticleCombo* locParticleCombo)
{
	//if is duplicate, returns the newly-created combo from the matching previous particle combo
	//otherwise, returns NULL
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	for(auto locNewComboPair : locNewParticleComboMap)
	{
		const DParticleCombo* locPreviousParticleCombo = locNewComboPair.first;
		if(locReaction == locPreviousParticleCombo->Get_Reaction()) //probably shouldn't be possible
			continue; //dreaction is the same: particle combos were different for a reason, even if the kinfit results are the same: keep both

		//see if steps are identical //may have an omega or phi resonance
		if(!locReaction->Check_AreStepsIdentical(locPreviousParticleCombo->Get_Reaction()))
			continue; //steps are not identical
		//see if particles are identical //particles may be in a different order
		if(!locParticleCombo->Check_AreMeasuredParticlesIdentical(locPreviousParticleCombo))
			continue; //particles are not identical

		return locNewComboPair.second;
	}

	return NULL;
}

void DParticleCombo_factory::Set_DecayingParticles(const DParticleCombo* locNewParticleCombo, const DParticleCombo* locOldParticleCombo, size_t locStepIndex, DParticleComboStep* locNewParticleComboStep, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults)
{
	DKinFitParticle* locKinFitParticle = Get_DecayingParticle(locOldParticleCombo, locStepIndex, locKinFitChain, locKinFitResults);
	if(locKinFitParticle == NULL) //not used in fit
	{
		locNewParticleComboStep->Set_InitialParticle(NULL);
		return; //no need to back-set NULL: was set to NULL by default
	}

	Particle_t locPID = PDGtoPType(locKinFitParticle->Get_PID());
	DKinematicData* locKinematicData_Position = Build_KinematicData(locPID, locKinFitParticle);
	DKinematicData* locKinematicData_Common = Build_KinematicData(locPID, locKinFitParticle);
	if(locKinFitParticle->Get_CommonVxParamIndex() >= 0)
		dKinFitUtils->Propagate_TrackInfoToCommonVertex(locKinematicData_Common, locKinFitParticle, &locKinFitResults->Get_VXi());

	bool locAtProdVertexFlag = locKinFitParticle->Get_VertexP4AtProductionVertex();
	DKinematicData* locKinematicData_InitState = locAtProdVertexFlag ? locKinematicData_Common : locKinematicData_Position;
	DKinematicData* locKinematicData_FinalState = locAtProdVertexFlag ? locKinematicData_Position : locKinematicData_Common;

	locNewParticleComboStep->Set_InitialParticle(locKinematicData_InitState);

	//now, back-set the particle at the other vertex
	int locFromStepIndex = locNewParticleComboStep->Get_InitialParticleDecayFromStepIndex();
	if((locStepIndex == 0) || (locFromStepIndex < 0))
		return; //no other place to set it

	DParticleComboStep* locParticleComboStep = const_cast<DParticleComboStep*>(locNewParticleCombo->Get_ParticleComboStep(locFromStepIndex));
	//find where it is the decaying particle
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Get_DecayStepIndex(loc_i) != int(locStepIndex))
			continue;
		locParticleComboStep->Set_FinalParticle(locKinematicData_FinalState, loc_i);
		break;
	}
}

DKinFitParticle* DParticleCombo_factory::Get_DecayingParticle(const DParticleCombo* locOldParticleCombo, size_t locComboStepIndex, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults)
{
	const DParticleComboStep* locParticleComboStep = locOldParticleCombo->Get_ParticleComboStep(locComboStepIndex);
	Particle_t locPID = locParticleComboStep->Get_InitialParticleID();
	if(!IsFixedMass(locPID))
		return NULL;

	//find which step in the DKinFitChain this combo step corresponds to
	for(size_t loc_i = 0; loc_i < locKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
	{
		const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(loc_i);

		//loop over init particles to get the decaying particle (if present)
		DKinFitParticle* locDecayingParticle = NULL;
		set<DKinFitParticle*> locInitialParticles = locKinFitChainStep->Get_InitialParticles();
		set<DKinFitParticle*>::iterator locParticleIterator = locInitialParticles.begin();
		for(; locParticleIterator != locInitialParticles.end(); ++locParticleIterator)
		{
			if((*locParticleIterator)->Get_KinFitParticleType() != d_DecayingParticle)
				continue; //not a decaying particle
			locDecayingParticle = *locParticleIterator;
			break;
		}
		if(locDecayingParticle == NULL)
			continue; //no decaying particles in this step

		if(PDGtoPType(locDecayingParticle->Get_PID()) != locPID)
			continue; //wrong PID

		//may still not be correct particle. compare decay products: if any of the combo step particles are in the kinfit step, this is it
			//if all step final particles are decaying, then dive down: through steps
			//if any step decay product at any step is located as any decay product at any step in the kinfit chain: then matches

		//get all measured products, then just pick the first one to search for
		deque<const DKinematicData*> locMeasuredParticles;
		locOldParticleCombo->Get_DecayChainParticles_Measured(locComboStepIndex, locMeasuredParticles);
		const DKinematicData* locMeasuredParticle = locMeasuredParticles[0];
		DKinFitParticle* locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locMeasuredParticle);
		if(locKinFitParticle == NULL) //null: neutral shower. Use shower object
		{
			const DNeutralShower* locNeutralShower = NULL;
			locMeasuredParticle->GetSingle(locNeutralShower);
			locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locNeutralShower);
		}

		if(!Search_ForParticleInDecay(locKinFitChain, loc_i, locKinFitParticle))
			continue;

		//Found!
		return locDecayingParticle;
	}

	return NULL;
}

bool DParticleCombo_factory::Search_ForParticleInDecay(const DKinFitChain* locKinFitChain, size_t locStepToSearch, DKinFitParticle* locParticleToFind)
{
	const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(locStepToSearch);
	set<DKinFitParticle*> locFinalParticles = locKinFitChainStep->Get_FinalParticles();
	if(locFinalParticles.find(locParticleToFind) != locFinalParticles.end())
		return true; //found it

	//else loop over final state, diving through decays
	set<DKinFitParticle*>::iterator locParticleIterator = locFinalParticles.begin();
	for(; locParticleIterator != locFinalParticles.end(); ++locParticleIterator)
	{
		if((*locParticleIterator)->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //not a decaying particle

		int locDecayStepIndex = locKinFitChain->Get_DecayStepIndex(*locParticleIterator);
		if(Search_ForParticleInDecay(locKinFitChain, locDecayStepIndex, locParticleToFind))
			return true; //found in in subsequent step
	}
	return false; //not found (yet)
}

const DChargedTrackHypothesis* DParticleCombo_factory::Get_ChargedHypothesis(const DParticleCombo* locParticleCombo, const vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses, const DKinematicData* locKinematicData_Measured) const
{
	for(size_t loc_l = 0; loc_l < locChargedTrackHypotheses.size(); ++loc_l)
	{
		if(!locChargedTrackHypotheses[loc_l]->IsAssociated(locParticleCombo))
			continue; // track created for a different particle combo
		if(locChargedTrackHypotheses[loc_l]->IsAssociated(locKinematicData_Measured))
			return locChargedTrackHypotheses[loc_l];
	}
	return dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData_Measured); //not used in fit: re-set the measured
}

const DNeutralParticleHypothesis* DParticleCombo_factory::Get_NeutralHypothesis(const DParticleCombo* locParticleCombo, const vector<const DNeutralParticleHypothesis*>& locNeutralParticleHypotheses, const DKinematicData* locKinematicData_Measured) const
{
	for(size_t loc_l = 0; loc_l < locNeutralParticleHypotheses.size(); ++loc_l)
	{
		if(!locNeutralParticleHypotheses[loc_l]->IsAssociated(locParticleCombo))
			continue; // track created for a different particle combo
		if(locNeutralParticleHypotheses[loc_l]->IsAssociated(locKinematicData_Measured))
			return locNeutralParticleHypotheses[loc_l];
	}
	return dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData_Measured); //not used in fit: re-set the measured
}

void DParticleCombo_factory::Set_SpacetimeVertex(const DParticleCombo* locNewParticleCombo, DParticleComboStep* locNewParticleComboStep, size_t locStepIndex, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain) const
{
	DKinFitType locKinFitType = locNewParticleCombo->Get_Reaction()->Get_KinFitType();
	if((locKinFitType == d_NoFit) || (locKinFitType == d_P4Fit))
		return; //neither vertex nor time was fit: no update to give

	//Position & Time
	const DKinematicData* locKinematicData = locNewParticleComboStep->Get_InitialParticle();
	if(locKinematicData != NULL)
	{
		locNewParticleComboStep->Set_Position(locKinematicData->position());

		bool locWasTimeKinFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));
		if(locWasTimeKinFitFlag)
			locNewParticleComboStep->Set_Time(locKinematicData->time());

		return;
	}

	//mass not fixed: if not initial step, get spacetime vertex from step where this particle was produced
	if(locStepIndex != 0)
	{
		//get spacetime vertex from step where this particle was produced
		int locDecayFromStepIndex = locNewParticleComboStep->Get_InitialParticleDecayFromStepIndex();
		const DParticleComboStep* locPreviousParticleComboStep = locNewParticleCombo->Get_ParticleComboStep(locDecayFromStepIndex);
		locNewParticleComboStep->Set_SpacetimeVertex(locPreviousParticleComboStep->Get_SpacetimeVertex());
		return;
	}

	//instead, get from common vertex of final state particles
	DKinFitParticle* locFinalKinFitParticle = *(locKinFitChain->Get_KinFitChainStep(0)->Get_FinalParticles().begin());

	//need the spacetime vertex at the production vertex of the particle grabbed
	TLorentzVector locSpacetimeVertex;
	if(locFinalKinFitParticle->Get_VertexP4AtProductionVertex()) //"position" is at production vertex
		locSpacetimeVertex = locFinalKinFitParticle->Get_SpacetimeVertex();
	else //"position" is at decay vertex
		locSpacetimeVertex = locFinalKinFitParticle->Get_CommonSpacetimeVertex(); //get production
	locNewParticleComboStep->Set_SpacetimeVertex(locSpacetimeVertex);
}

void DParticleCombo_factory::Reset_Data(void)
{
	//delete objects that this factory created (since the NOT_OBJECT_OWNER flag is set)
		//if there are no kinfit results then they were simply copied from the other factory
	for(size_t loc_i = 0; loc_i < dCreatedParticleCombos.size(); ++loc_i)
		delete dCreatedParticleCombos[loc_i];
	_data.clear();
	dCreatedParticleCombos.clear();
}

DKinematicData* DParticleCombo_factory::Build_KinematicData(Particle_t locPID, DKinFitParticle* locKinFitParticle)
{
	DKinematicData* locKinematicData = Get_KinematicDataResource();
	locKinematicData->setPID(locPID);
	locKinematicData->setCharge(ParticleCharge(locPID));
	if(locPID != Unknown)
		locKinematicData->setMass(ParticleMass(locPID));
	else
		locKinematicData->setMass(locKinFitParticle->Get_Mass());
	locKinematicData->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locKinematicData->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locKinematicData->setTime(locKinFitParticle->Get_Time());
	if(locKinFitParticle->Get_CovarianceMatrix() != NULL)
		locKinematicData->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());
	locKinematicData->setPathLength(locKinFitParticle->Get_PathLength(), locKinFitParticle->Get_PathLengthUncertainty());

	return locKinematicData;
}

void DParticleCombo_factory::Reset_Pools(void)
{
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dParticleComboStepPool_All.size() > MAX_DParticleComboStepPoolSize)
	{
		for(size_t loc_i = MAX_DParticleComboStepPoolSize; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
			delete dParticleComboStepPool_All[loc_i];
		dParticleComboStepPool_All.resize(MAX_DParticleComboStepPoolSize);
	}
	dParticleComboStepPool_Available = dParticleComboStepPool_All;

	if(dKinematicDataPool_All.size() > MAX_DKinematicDataPoolSize)
	{
		for(size_t loc_i = MAX_DKinematicDataPoolSize; loc_i < dKinematicDataPool_All.size(); ++loc_i)
			delete dKinematicDataPool_All[loc_i];
		dKinematicDataPool_All.resize(MAX_DKinematicDataPoolSize);
	}
	dKinematicDataPool_Available = dKinematicDataPool_All;
}

DKinematicData* DParticleCombo_factory::Get_KinematicDataResource(void)
{
	DKinematicData* locKinematicData;
	if(dKinematicDataPool_Available.empty())
	{
		locKinematicData = new DKinematicData;
		dKinematicDataPool_All.push_back(locKinematicData);
	}
	else
	{
		locKinematicData = dKinematicDataPool_Available.back();
		locKinematicData->Reset();
		dKinematicDataPool_Available.pop_back();
	}
	return locKinematicData;
}

DParticleComboStep* DParticleCombo_factory::Get_ParticleComboStepResource(void)
{
	DParticleComboStep* locParticleComboStep;
	if(dParticleComboStepPool_Available.empty())
	{
		locParticleComboStep = new DParticleComboStep;
		dParticleComboStepPool_All.push_back(locParticleComboStep);
	}
	else
	{
		locParticleComboStep = dParticleComboStepPool_Available.back();
		locParticleComboStep->Reset();
		dParticleComboStepPool_Available.pop_back();
	}
	return locParticleComboStep;
}

//------------------
// erun
//------------------
jerror_t DParticleCombo_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleCombo_factory::fini(void)
{
	Reset_Data();

	for(size_t loc_i = 0; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
		delete dParticleComboStepPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinematicDataPool_All.size(); ++loc_i)
		delete dKinematicDataPool_All[loc_i];

	return NOERROR;
}
