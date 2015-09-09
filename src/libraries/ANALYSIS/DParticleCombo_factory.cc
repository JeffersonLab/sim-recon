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
jerror_t DParticleCombo_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	const DMagneticFieldMap* locMagneticFieldMap = locApplication->GetBfield(runnumber);

	double locTargetZCenter = 65.0;
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	locGeometry->GetTargetZ(locTargetZCenter);

	double locBx, locBy, locBz;
	locMagneticFieldMap->GetField(0.0, 0.0, locTargetZCenter, locBx, locBy, locBz);
	TVector3 locBField(locBx, locBy, locBz);
	if(locBField.Mag() > 0.0)
		dKinFitter.Set_BField(locMagneticFieldMap);

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
jerror_t DParticleCombo_factory::evnt(JEventLoop* locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DParticleCombo_factory::evnt()");
#endif

	Reset_Data();
	Reset_Pools();

	vector<const DKinFitResults*> locKinFitResultsVector;
	locEventLoop->Get(locKinFitResultsVector);

 	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector, "PreKinFit");

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
		//all failed! save 'em and bail
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

	DKinematicData* locKinematicData;
	const DKinFitParticle* locKinFitParticle;
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

			locParticleCombos_FailedKinFit.erase(locParticleCombo); //kinfit successful, don't copy pointer later

			DParticleCombo* locNewParticleCombo = new DParticleCombo();
			dCreatedParticleCombos.push_back(locNewParticleCombo);
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
				locNewParticleComboStep->Set_MeasuredParticleComboStep(locParticleComboStep);
				locNewParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboStep->Get_ParticleComboBlueprintStep());
				locNewParticleComboStep->Set_SpacetimeVertex(locParticleComboStep->Get_SpacetimeVertex()); //overridden if kinematic fit
				bool locWasVertexKinFitFlag = ((locKinFitType != d_NoFit) && (locKinFitType != d_P4Fit));
				bool locWasTimeKinFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));

				//INITIAL PARTICLE & SPACETIME VERTEX
				locPID = locParticleComboStep->Get_InitialParticleID();
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
					const DKinematicData* locKinematicData_Measured = locParticleComboStep->Get_FinalParticle_Measured(loc_k);
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
							vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses_Associated;
							locNeutralParticleHypotheses[loc_l]->Get(locNeutralParticleHypotheses_Associated);
							//loop over associated: default tag, "Combo" tag, etc.
							bool locMatchFlag = false;
							for(size_t loc_m = 0; loc_m < locNeutralParticleHypotheses_Associated.size(); ++loc_m)
							{
								if(locNeutralParticleHypotheses_Associated[loc_m] != locKinematicData_Measured)
									continue; //wrong neutral hypothesis
								locNeutralParticleHypotheses[loc_l]->Get(locParticleCombos_Associated);
								for(size_t loc_n = 0; loc_n < locParticleCombos_Associated.size(); ++loc_n)
								{
									if(locParticleCombos_Associated[loc_n] != locParticleCombo)
										continue;
									locMatchFlag = true;
									break;
								}
								if(!locMatchFlag)
									continue; // track created for a different particle combo
								locNewParticleComboStep->Add_FinalParticle(locNeutralParticleHypotheses[loc_l]);
								break;
							}
							if(locMatchFlag)
								break;
						}
					}
					else //charged
					{
						for(size_t loc_l = 0; loc_l < locChargedTrackHypotheses.size(); ++loc_l)
						{
							vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses_Associated;
							locChargedTrackHypotheses[loc_l]->Get(locChargedTrackHypotheses_Associated);
							//loop over associated: default tag, "Combo" tag, etc.
							bool locMatchFlag = false;
							for(size_t loc_m = 0; loc_m < locChargedTrackHypotheses_Associated.size(); ++loc_m)
							{
								if(locChargedTrackHypotheses_Associated[loc_m] != locKinematicData_Measured)
									continue; //wrong track hypothesis
								locChargedTrackHypotheses[loc_l]->Get(locParticleCombos_Associated);
								for(size_t loc_n = 0; loc_n < locParticleCombos_Associated.size(); ++loc_n)
								{
									if(locParticleCombos_Associated[loc_n] != locParticleCombo)
										continue;
									locMatchFlag = true;
									break;
								}
								if(!locMatchFlag)
									continue; // track created for a different particle combo
								locNewParticleComboStep->Add_FinalParticle(locChargedTrackHypotheses[loc_l]);
								break;
							}
							if(locMatchFlag)
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

	//directly save all combos for which the kinfits failed
	set<const DParticleCombo*>::iterator locIterator = locParticleCombos_FailedKinFit.begin();
	for(; locIterator != locParticleCombos_FailedKinFit.end(); ++locIterator)
		_data.push_back(const_cast<DParticleCombo*>(*locIterator));

	return NOERROR;
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
	for(size_t loc_i = 0; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
		delete dParticleComboStepPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinematicDataPool_All.size(); ++loc_i)
		delete dKinematicDataPool_All[loc_i];

	return NOERROR;
}

