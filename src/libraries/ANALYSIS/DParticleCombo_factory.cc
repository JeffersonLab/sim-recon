#include "DParticleCombo_factory.h"

//------------------
// init
//------------------
jerror_t DParticleCombo_factory::init(void)
{
	MAX_DParticleComboStepPoolSize = 40;
	MAX_DKinematicDataPoolSize = 40;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleCombo_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	dKinFitter.Set_BField(locApplication->GetBfield());

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

	MAX_DParticleComboStepPoolSize = 2000*locNumReactions;
	MAX_DKinematicDataPoolSize = 1000*locNumReactions;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleCombo_factory::evnt(JEventLoop* locEventLoop, int eventnumber)
{
	Reset_Pools();

	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

 	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector, "PreKinFit");

	deque<const DParticleCombo*> locParticleCombos_PreKinFit, locSurvivedParticleCombos;
	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		locAnalysisResultsVector[loc_i]->Get_PassedParticleCombos(locSurvivedParticleCombos);
		locParticleCombos_PreKinFit.insert(locParticleCombos_PreKinFit.end(), locSurvivedParticleCombos.begin(), locSurvivedParticleCombos.end());
	}

	if(locKinFitResultsVector.empty())
	{
		//kinfit not requested (or all kinfits failed), just clone the original objects
		for(size_t loc_i = 0; loc_i < locParticleCombos_PreKinFit.size(); ++loc_i)
			_data.push_back(Clone_ParticleCombo(locParticleCombos_PreKinFit[loc_i]));
		return NOERROR;
	}

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses, "KinFit");

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses, "KinFit");

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons, "KinFit");

	DKinematicData* locKinematicData;
	const DKinematicData* locConstKinematicData;
	const DKinFitParticle* locKinFitParticle;
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses_Associated;
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses_Associated;
	vector<const DParticleComboBlueprint*> locParticleComboBlueprints_Associated;
	vector<const DParticleCombo*> locParticleCombos_Associated;
	Particle_t locPID;
	map<pair<Particle_t, deque<const DKinematicData*> >, const DKinFitParticle*> locDecayingParticles;

	map<const DParticleCombo*, DParticleCombo*> locKinFitParticleComboMap;
	for(size_t loc_i = 0; loc_i < locKinFitResultsVector.size(); ++loc_i)
	{
		DKinFitType locKinFitType = locKinFitResultsVector[loc_i]->Get_KinFitType();

		set<const DParticleCombo*> locParticleCombos;
		locKinFitResultsVector[loc_i]->Get_ParticleCombos(locParticleCombos); //each result may be for several combos (duplicate results)

		set<const DParticleCombo*>::iterator locComboIterator = locParticleCombos.begin();
		for(; locComboIterator != locParticleCombos.end(); ++locComboIterator)
		{
			const DParticleCombo* locParticleCombo = *locComboIterator;
			const DReaction* locReaction = locParticleCombo->Get_Reaction();

			DParticleCombo* locNewParticleCombo = new DParticleCombo();
			locNewParticleCombo->Set_Reaction(locParticleCombo->Get_Reaction());
			locNewParticleCombo->Set_KinFitResults(locKinFitResultsVector[loc_i]);
			locNewParticleCombo->Set_EventRFBunch(locParticleCombo->Get_EventRFBunch());

			locParticleCombo->GetT(locParticleComboBlueprints_Associated);
			locNewParticleCombo->AddAssociatedObject(locParticleComboBlueprints_Associated[0]);
			locNewParticleCombo->AddAssociatedObject(locParticleCombo);

			//search to see if the results were the same as a previous combo. if so, copy them:
				//note that it is possible for two combos with different steps to have the same kinfit results (e.g. one has an omega or phi, and the other doesn't)
			bool locMatchFlag = false;
			set<const DParticleCombo*>::iterator locPreviousComboIterator = locParticleCombos.begin();
			for(; locPreviousComboIterator != locComboIterator; ++locPreviousComboIterator)
			{
				const DParticleCombo* locPreviousParticleCombo = *locPreviousComboIterator;
				if(locParticleCombo->Get_Reaction() == locPreviousParticleCombo->Get_Reaction()) //probably shouldn't be possible
					continue; //dreaction is the same: particle combos were different for a reason, even if the kinfit results are the same (e.g. ???): keep both

				//see if steps are identical //may have an omega or phi resonance
				if(!locReaction->Check_AreStepsIdentical(locPreviousParticleCombo->Get_Reaction()))
					continue; //steps are not identical
				//see if particles are identical //particles may be in a different order
				if(!locParticleCombo->Check_AreMeasuredParticlesIdentical(locPreviousParticleCombo))
					continue; //particles are not identical

				//everything is identical: copy results
				const DParticleCombo* locPreviousNewParticleCombo = locKinFitParticleComboMap[locPreviousParticleCombo];
				for(size_t loc_k = 0; loc_k < locPreviousNewParticleCombo->Get_NumParticleComboSteps(); ++loc_k)
					locNewParticleCombo->Add_ParticleComboStep(locPreviousNewParticleCombo->Get_ParticleComboStep(loc_k));
				locMatchFlag = true;
				break;
			}
			if(locMatchFlag)
			{
				_data.push_back(locNewParticleCombo);
				continue;
			}

			//unique combo: build steps
			locKinFitResultsVector[loc_i]->Get_DecayingParticles(locDecayingParticles);
			for(size_t loc_j = 0; loc_j < locParticleCombo->Get_NumParticleComboSteps(); ++loc_j)
			{
				const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_j);
				DParticleComboStep* locNewParticleComboStep = Get_ParticleComboStepResource();
				locNewParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboStep->Get_ParticleComboBlueprintStep());
				locNewParticleComboStep->Set_SpacetimeVertex(locParticleComboStep->Get_SpacetimeVertex()); //overridden if kinematic fit
				bool locWasVertexKinFitFlag = ((locKinFitType != d_NoFit) && (locKinFitType != d_P4Fit));
				bool locWasTimeKinFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));

				//INITIAL PARTICLE & SPACETIME VERTEX
				locPID = locParticleComboStep->Get_InitialParticleID();
				locNewParticleComboStep->Set_InitialParticle_Measured(locParticleComboStep->Get_InitialParticle_Measured());
				if(locParticleComboStep->Is_InitialParticleDetected()) //set beam photon
				{
					for(size_t loc_k = 0; loc_k < locBeamPhotons.size(); ++loc_k)
					{
						locBeamPhotons[loc_k]->GetT(locParticleCombos_Associated);
						bool locMatchFlag = false;
						for(size_t loc_m = 0; loc_m < locParticleCombos_Associated.size(); ++loc_m)
						{
							if(locParticleCombos_Associated[loc_m] != locParticleCombo)
								continue;
							locMatchFlag = true;
							break;
						}
						if(!locMatchFlag)
							continue; // beam photon created for a different particle combo
						locNewParticleComboStep->Set_InitialParticle(locBeamPhotons[loc_k]);
						if(locWasVertexKinFitFlag)
							locNewParticleComboStep->Set_Position(locBeamPhotons[loc_k]->position());
						if(locWasTimeKinFitFlag)
							locNewParticleComboStep->Set_Time(locBeamPhotons[loc_k]->time());
						break;
					}
				}
				else //decaying particle!
				{
					deque<const DKinematicData*> locDecayProducts;
					locParticleCombo->Get_DecayChainParticles_Measured(loc_j, locDecayProducts);
					pair<Particle_t, deque<const DKinematicData*> > locDecayInfoPair(locPID, locDecayProducts);
					locKinFitParticle = locDecayingParticles[locDecayInfoPair];
					if(locKinFitParticle != NULL)
					{
						locKinematicData = Build_KinematicData(locPID, locKinFitParticle);
						//if true: propagate the track info from the production vertex to the decay vertex
						if(locKinFitParticle->Get_DecayingParticleAtProductionVertexFlag() && (locKinFitParticle->Get_NumVertexFits() == 2))
							dKinFitter.Propagate_TrackInfoToCommonVertex(locKinematicData, locKinFitParticle, locKinFitResultsVector[loc_i]->Get_VXi());
						if(locWasVertexKinFitFlag)
							locNewParticleComboStep->Set_Position(locKinematicData->position());
						if(locWasTimeKinFitFlag)
							locNewParticleComboStep->Set_Time(locKinematicData->time());
					}
					else
					{
						locKinematicData = NULL;
						//if initial step, or vertex not fit, get from pre-kinfit
						//resonance: get spacetime vertex from step where this particle was produced
						int locInitialParticleDecayFromStepIndex = locParticleComboStep->Get_InitialParticleDecayFromStepIndex();
						const DParticleComboStep* locPreviousParticleComboStep = locParticleCombo->Get_ParticleComboStep(locInitialParticleDecayFromStepIndex);
						locNewParticleComboStep->Set_SpacetimeVertex(locPreviousParticleComboStep->Get_SpacetimeVertex());
					}
					locNewParticleComboStep->Set_InitialParticle(locKinematicData);
				}

				//TARGET PARTICLE //set position and time??
				locPID = locParticleComboStep->Get_InitialParticleID();
				if(locParticleComboStep->Is_TargetPresent())
					locNewParticleComboStep->Set_TargetParticle(locParticleComboStep->Get_TargetParticle());

				//FINAL PARTICLES
				for(size_t loc_k = 0; loc_k < locParticleComboStep->Get_NumFinalParticles(); ++loc_k)
				{
					locConstKinematicData = locParticleComboStep->Get_FinalParticle_Measured(loc_k);
					locNewParticleComboStep->Add_FinalParticle_Measured(locConstKinematicData);
					locPID = locParticleComboStep->Get_FinalParticleID(loc_k);
					if(locParticleComboStep->Is_FinalParticleMissing(loc_k)) //missing!
					{
						locKinFitParticle = locKinFitResultsVector[loc_i]->Get_MissingParticle();;
						locKinematicData = (locKinFitParticle != NULL) ? Build_KinematicData(locPID, locKinFitParticle) : NULL;
						locNewParticleComboStep->Add_FinalParticle(locKinematicData);
					}
					else if(locParticleComboStep->Is_FinalParticleDecaying(loc_k)) //decaying!
					{
						deque<const DKinematicData*> locDecayProducts;
						locParticleCombo->Get_DecayChainParticles_Measured(locParticleComboStep->Get_DecayStepIndex(loc_k), locDecayProducts);
						pair<Particle_t, deque<const DKinematicData*> > locDecayInfoPair(locPID, locDecayProducts);
						locKinFitParticle = locDecayingParticles[locDecayInfoPair];
						if(locKinFitParticle != NULL)
						{
							locKinematicData = Build_KinematicData(locPID, locKinFitParticle);
							//if true: propagate the track info from the decay vertex to the production vertex
							if(!locKinFitParticle->Get_DecayingParticleAtProductionVertexFlag() && (locKinFitParticle->Get_NumVertexFits() == 2))
								dKinFitter.Propagate_TrackInfoToCommonVertex(locKinematicData, locKinFitParticle, locKinFitResultsVector[loc_i]->Get_VXi());
						}
						else
							locKinematicData = NULL;
						locNewParticleComboStep->Add_FinalParticle(locKinematicData);
					}
					else if(locParticleComboStep->Is_FinalParticleNeutral(loc_k)) //neutral
					{
						for(size_t loc_l = 0; loc_l < locNeutralParticleHypotheses.size(); ++loc_l)
						{
							locNeutralParticleHypotheses[loc_l]->GetT(locNeutralParticleHypotheses_Associated);
							if(locNeutralParticleHypotheses_Associated[0] != locConstKinematicData)
								continue; //wrong track hypothesis
							locNeutralParticleHypotheses[loc_l]->GetT(locParticleCombos_Associated);
							bool locMatchFlag = false;
							for(size_t loc_m = 0; loc_m < locParticleCombos_Associated.size(); ++loc_m)
							{
								if(locParticleCombos_Associated[loc_m] != locParticleCombo)
									continue;
								locMatchFlag = true;
								break;
							}
							if(!locMatchFlag)
								continue; // track created for a different particle combo
							locNewParticleComboStep->Add_FinalParticle(locNeutralParticleHypotheses[loc_l]);
							break;
						}
					}
					else //charged
					{
						for(size_t loc_l = 0; loc_l < locChargedTrackHypotheses.size(); ++loc_l)
						{
							locChargedTrackHypotheses[loc_l]->GetT(locChargedTrackHypotheses_Associated);
							if(locChargedTrackHypotheses_Associated[0] != locConstKinematicData)
								continue; //wrong track hypothesis
							locChargedTrackHypotheses[loc_l]->GetT(locParticleCombos_Associated);
							bool locMatchFlag = false;
							for(size_t loc_m = 0; loc_m < locParticleCombos_Associated.size(); ++loc_m)
							{
								if(locParticleCombos_Associated[loc_m] != locParticleCombo)
									continue;
								locMatchFlag = true;
								break;
							}
							if(!locMatchFlag)
								continue; // track created for a different particle combo
							locNewParticleComboStep->Add_FinalParticle(locChargedTrackHypotheses[loc_l]);
							break;
						}
					}
				}
				locNewParticleCombo->Add_ParticleComboStep(locNewParticleComboStep);
			}

			locKinFitParticleComboMap[locParticleCombo] = locNewParticleCombo;
			_data.push_back(locNewParticleCombo);
		}
	}

	//clone all combos for which the kinfits failed or didn't occur
	for(size_t loc_i = 0; loc_i < _data.size(); ++loc_i) //first remove successfully kinfit combos from the locParticleCombos_PreKinFit list
	{
		_data[loc_i]->GetT(locParticleCombos_Associated);
		for(deque<const DParticleCombo*>::iterator locIterator = locParticleCombos_PreKinFit.begin(); locIterator != locParticleCombos_PreKinFit.end(); ++locIterator)
		{
			if((*locIterator) == locParticleCombos_Associated[0])
			{
				locParticleCombos_PreKinFit.erase(locIterator); //kinfit did not fail
				break;
			}
		}
	}
	for(size_t loc_i = 0; loc_i < locParticleCombos_PreKinFit.size(); ++loc_i) //kinfit failed for these, clone and save new ones
	{
		DParticleCombo* locNewParticleCombo = Clone_ParticleCombo(locParticleCombos_PreKinFit[loc_i]);
		_data.push_back(locNewParticleCombo);
	}

	return NOERROR;
}

DParticleCombo* DParticleCombo_factory::Clone_ParticleCombo(const DParticleCombo* locParticleCombo)
{
	DParticleCombo* locNewParticleCombo = new DParticleCombo(*locParticleCombo);
	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locParticleCombo->GetT(locParticleComboBlueprints);
	locNewParticleCombo->AddAssociatedObject(locParticleComboBlueprints[0]);
	return locNewParticleCombo;
}

DKinematicData* DParticleCombo_factory::Build_KinematicData(Particle_t locPID, const DKinFitParticle* locKinFitParticle)
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
		locKinematicData->setErrorMatrix(*(locKinFitParticle->Get_CovarianceMatrix()));
	locKinematicData->setPathLength(locKinFitParticle->Get_PathLength(), locKinFitParticle->Get_PathLengthUncertainty());

	return locKinematicData;
}

void DParticleCombo_factory::Reset_Pools(void)
{
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dParticleComboStepPool_All.size() > MAX_DParticleComboStepPoolSize){
		for(size_t loc_i = MAX_DParticleComboStepPoolSize; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
			delete dParticleComboStepPool_All[loc_i];
		dParticleComboStepPool_All.resize(MAX_DParticleComboStepPoolSize);
	}
	dParticleComboStepPool_Available = dParticleComboStepPool_All;

	if(dKinematicDataPool_All.size() > MAX_DKinematicDataPoolSize){
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
		Reset_KinematicData(locKinematicData);
		locKinematicData->ClearAssociatedObjects();
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

void DParticleCombo_factory::Reset_KinematicData(DKinematicData* locKinematicData)
{
	locKinematicData->setPID(Unknown);
	locKinematicData->setMassFixed();
	locKinematicData->setCharge(0);
	locKinematicData->setMass(0.0);

	locKinematicData->setMomentum(DVector3());
	locKinematicData->setPosition(DVector3());
	locKinematicData->setTime(0.0);

	locKinematicData->setdEdx(0.0);
	locKinematicData->setPathLength(0.0, 0.0);
	locKinematicData->setTrackingStateVector(0.0, 0.0, 0.0, 0.0, 0.0);

	locKinematicData->setT0(0.0, 0.0, SYS_NULL);
	locKinematicData->setT1(0.0, 0.0, SYS_NULL);

	locKinematicData->clearErrorMatrix();
	locKinematicData->clearTrackingErrorMatrix();
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
	return NOERROR;
}

