#include "DKinFitResults_factory.h"

//------------------
// init
//------------------
jerror_t DKinFitResults_factory::init(void)
{
	dDebugLevel = 0;
	dKinFitDebugLevel = 0;
	dLinkVerticesFlag = true;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DKinFitResults_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	if(locAnalysisUtilitiesVector.empty())
	{
		_DBG_<<"Unable to get a DAnalysisUtilities object!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	const DMagneticFieldMap* locMagneticFieldMap = locApplication->GetBfield();
	dKinFitter.Set_BField(locMagneticFieldMap);

	gPARMS->SetDefaultParameter("KINFIT:KINFITDEBUGLEVEL", dKinFitDebugLevel);
	gPARMS->SetDefaultParameter("KINFIT:DEBUGLEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);

	dKinFitter.Set_DebugLevel(dKinFitDebugLevel);
	dKinFitter.Set_LinkVerticesFlag(dLinkVerticesFlag);

	dTargetZCenter = 65.0;
	// Get Target parameters from XML
	DGeometry* locGeometry = locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	if(locGeometry != NULL)
		locGeometry->GetTargetZ(dTargetZCenter);

	//test
	double locBx, locBy, locBz;
	locMagneticFieldMap->GetField(0.0, 0.0, dTargetZCenter, locBx, locBy, locBz);
	TVector3 locBField(locBx, locBy, locBz);
	if(!(locBField.Mag() > 0.0))
		cout << "WARNING: MAGNETIC FIELD IS ZERO AT THE TARGET CENTER. YOU SURE THIS IS OK???" << endl;

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

	//set pool sizes
	unsigned int locExpectedNumCombos = 1000;
	dKinFitter.Set_MaxKinFitParticlePoolSize(locNumReactions*locExpectedNumCombos*6);
	dKinFitter.Set_MaxKinFitConstraintVertexPoolSize(locNumReactions*locExpectedNumCombos*3);
	dKinFitter.Set_MaxKinFitConstraintSpacetimePoolSize(locNumReactions*locExpectedNumCombos*3);
	dKinFitter.Set_MaxKinFitConstraintP4PoolSize(locNumReactions*locExpectedNumCombos*3);
	dKinFitter.Set_MaxMatrixDSymPoolSize(locNumReactions*locExpectedNumCombos*6);

	return NOERROR;
}

void DKinFitResults_factory::Reset_NewEvent(void)
{
	dKinFitter.Reset_NewEvent();
}

//------------------
// evnt
//------------------
jerror_t DKinFitResults_factory::evnt(JEventLoop* locEventLoop, int eventnumber)
{
	dPreviouslyFailedFits.clear();

	//perform all of the analysis steps that don't need the kinematic fit results (saves time by reducing #kinfits)
 	vector<const DAnalysisResults*> locAnalysisResultsVector;
	locEventLoop->Get(locAnalysisResultsVector, "PreKinFit");

	//get all of the ParticleCombos that survive the cuts
	deque<const DParticleCombo*> locParticleCombos, locSurvivedParticleCombos;
	for(size_t loc_i = 0; loc_i < locAnalysisResultsVector.size(); ++loc_i)
	{
		locAnalysisResultsVector[loc_i]->Get_PassedParticleCombos(locSurvivedParticleCombos);
		locParticleCombos.insert(locParticleCombos.end(), locSurvivedParticleCombos.begin(), locSurvivedParticleCombos.end());
	}

	dKinFitter.Reset_NewEvent();
	for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
	{
		const DParticleCombo* locParticleCombo = locParticleCombos[loc_i];
		const DReaction* locReaction = locParticleCombo->Get_Reaction();
		if(locReaction->Get_KinFitType() == d_NoFit)
			continue; //don't do any kinematic fits!

		if(dDebugLevel > 0)
			cout << "Create kinematic fit constraints for event: " << eventnumber << endl;

		map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles;
		deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > > locSortedConstraints;
		deque<DKinFitConstraint*> locOriginalConstraints;
		if(!Create_KinFitConstraints(locParticleCombo, locDecayingKinFitParticles, locOriginalConstraints, locSortedConstraints))
			continue; //sort-constraints failed: invalid! cannot setup kinfit

		//ok to just grab the rf bunch from the first combo: if they are different between combos, then they weren't needed anyway
		const DEventRFBunch* locEventRFBunch = locParticleCombo->Get_EventRFBunch();

		//check if previous combos would have resulted in duplicate results: if so, just copy the result (or abort fit if it failed last time)
			//this could happen if you want to perform the exact same kinematic fit, but want to perform different DAnalysisActions on them
		if(Handle_IfKinFitResultsWillBeIdentical(locParticleCombo, locOriginalConstraints, locEventRFBunch, locDecayingKinFitParticles))
			continue; //identical results

		//setup kinfit
		DKinFitType locKinFitType = locParticleCombo->Get_Reaction()->Get_KinFitType();
		if(!Setup_KinFit(locKinFitType, locOriginalConstraints, locEventRFBunch, locSortedConstraints))
		{
			if(dDebugLevel > 0)
				cout << "Kinematic Fit Setup Failed." << endl;
			deque<const DKinFitConstraint*> locConstOriginalConstraints(locOriginalConstraints.begin(), locOriginalConstraints.end());
			dPreviouslyFailedFits.push_back(DPreviousFitInfo(locEventRFBunch, locConstOriginalConstraints, locDecayingKinFitParticles));
			continue;
		}
		if(dDebugLevel > 0)
			cout << "Perform Primary Kinematic Fit" << endl;
		if(dKinFitter.Fit_Reaction()) //if fit fails: no kinfit results, BUT will still generate new DParticleCombo (using old info though!)
		{
			if(dDebugLevel > 0)
				cout << "Kinematic Fit Converged for Event " << eventnumber << ", confidence = " << dKinFitter.Get_ConfidenceLevel() << endl;
			Build_KinFitResults(locParticleCombo, locDecayingKinFitParticles, locOriginalConstraints);
		}
		else
		{
			deque<const DKinFitConstraint*> locConstOriginalConstraints(locOriginalConstraints.begin(), locOriginalConstraints.end());
			dPreviouslyFailedFits.push_back(DPreviousFitInfo(locEventRFBunch, locConstOriginalConstraints, locDecayingKinFitParticles));
		}
	}

	return NOERROR;
}

bool DKinFitResults_factory::Handle_IfKinFitResultsWillBeIdentical(const DParticleCombo* locParticleCombo, deque<DKinFitConstraint*> locConstraints_ToCheck, const DEventRFBunch* locRFBunch_ToCheck, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_ToCheck)
{
	//first check previously-passed results
	for(size_t loc_j = 0; loc_j < _data.size(); ++loc_j)
	{
		deque<const DKinFitConstraint*> locPreviousOriginalConstraints;
		_data[loc_j]->Get_OriginalKinFitConstraints(locPreviousOriginalConstraints);

		set<const DParticleCombo*> locPreviousParticleCombos;
		_data[loc_j]->Get_ParticleCombos(locPreviousParticleCombos);
		const DEventRFBunch* locPreviousEventRFBunch = (*locPreviousParticleCombos.begin())->Get_EventRFBunch();

		map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_CheckAgainst;
		_data[loc_j]->Get_InputDecayingParticleInfo(locDecayingKinFitParticles_CheckAgainst);

		if(!Check_IfKinFitResultsWillBeIdentical(locDecayingKinFitParticles_ToCheck, locDecayingKinFitParticles_CheckAgainst))
			continue;
		if(!Check_IfKinFitResultsWillBeIdentical(locConstraints_ToCheck, locPreviousOriginalConstraints, locRFBunch_ToCheck, locPreviousEventRFBunch))
			continue;
		if(dDebugLevel > 0)
			cout << "kinfit results will be the same: register combo & skip fit" << endl;

		//kinfit results will be the same: clone and save
		_data[loc_j]->Add_ParticleCombo(locParticleCombo);
		return true;
	}

	//now check previously-failed fits
	for(size_t loc_j = 0; loc_j < dPreviouslyFailedFits.size(); ++loc_j)
	{
		if(!Check_IfKinFitResultsWillBeIdentical(locDecayingKinFitParticles_ToCheck, dPreviouslyFailedFits[loc_j].dDecayingParticles))
			continue;
		if(!Check_IfKinFitResultsWillBeIdentical(locConstraints_ToCheck, dPreviouslyFailedFits[loc_j].dOriginalConstraints, locRFBunch_ToCheck, dPreviouslyFailedFits[loc_j].dEventRFBunch))
			continue;
		if(dDebugLevel > 0)
			cout << "kinfit results will be the same & will fail: skip fit" << endl;
		return true;
	}
	return false;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_ToCheck, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles_CheckAgainst)
{
	if(locDecayingKinFitParticles_ToCheck.size() != locDecayingKinFitParticles_CheckAgainst.size())
		return false;

	map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >::iterator locIterator_ToCheck = locDecayingKinFitParticles_ToCheck.begin();
	for(; locIterator_ToCheck != locDecayingKinFitParticles_ToCheck.end(); ++locIterator_ToCheck)
	{
		map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >::iterator locIterator_CheckAgainst = locDecayingKinFitParticles_CheckAgainst.begin();
		bool locMatchFoundFlag = false;
		for(; locIterator_CheckAgainst != locDecayingKinFitParticles_CheckAgainst.end(); ++locIterator_CheckAgainst)
		{
			if(locIterator_ToCheck->second.first != locIterator_CheckAgainst->second.first)
				continue;
			deque<const DKinematicData*> locParticles_ToCheck = locIterator_ToCheck->second.second;
			deque<const DKinematicData*> locParticles_CheckAgainst = locIterator_CheckAgainst->second.second;
			if(locParticles_ToCheck.size() != locParticles_CheckAgainst.size())
				continue;

			//compare particles //note particles can be in a different order!
			bool locAllParticlesMatchFlag = true;
			for(size_t loc_i = 0; loc_i < locParticles_ToCheck.size(); ++loc_i)
			{
				bool locParticleMatchFoundFlag = false;
				deque<const DKinematicData*>::iterator locDequeIterator = locParticles_CheckAgainst.begin();
				for(; locDequeIterator != locParticles_CheckAgainst.end(); ++locDequeIterator)
				{
					if(locParticles_ToCheck[loc_i] != (*locDequeIterator))
						continue;
					locParticles_CheckAgainst.erase(locDequeIterator);
					locParticleMatchFoundFlag = true;
					break;
				}
				if(!locParticleMatchFoundFlag)
				{
					locAllParticlesMatchFlag = false;
					break;
				}
			}
			if(!locAllParticlesMatchFlag)
				continue;

			locMatchFoundFlag = true;
			locDecayingKinFitParticles_CheckAgainst.erase(locIterator_CheckAgainst);
			break;
		}
		if(!locMatchFoundFlag)
			return false;
	}
	return true;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(deque<DKinFitConstraint*> locConstraints_ToCheck, deque<const DKinFitConstraint*> locConstraints_CheckAgainst, const DEventRFBunch* locRFBunch_ToCheck, const DEventRFBunch* locRFBunch_CheckAgainst)
{
	if(locConstraints_ToCheck.size() != locConstraints_CheckAgainst.size())
		return false;

	for(size_t loc_i = 0; loc_i < locConstraints_ToCheck.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint_ToCheck = dynamic_cast<DKinFitConstraint_P4*>(locConstraints_ToCheck[loc_i]);
		DKinFitConstraint_Vertex* locVertexConstraint_ToCheck = dynamic_cast<DKinFitConstraint_Vertex*>(locConstraints_ToCheck[loc_i]);
		DKinFitConstraint_Spacetime* locSpacetimeConstraint_ToCheck = dynamic_cast<DKinFitConstraint_Spacetime*>(locConstraints_ToCheck[loc_i]);

		bool locMatchFoundFlag = false;
		deque<const DKinFitConstraint*>::iterator locIterator = locConstraints_CheckAgainst.begin();
		for(; locIterator != locConstraints_CheckAgainst.end(); ++locIterator)
		{
			if(locP4Constraint_ToCheck != NULL)
			{
				const DKinFitConstraint_P4* locP4Constraint_CheckAgainst = dynamic_cast<const DKinFitConstraint_P4*>(*locIterator);
				if(locP4Constraint_CheckAgainst == NULL)
					continue;
				if(!Check_IfKinFitResultsWillBeIdentical(locP4Constraint_ToCheck, locP4Constraint_CheckAgainst))
					continue;
				locMatchFoundFlag = true;
				locConstraints_CheckAgainst.erase(locIterator);
				break;
			}
			if(locVertexConstraint_ToCheck != NULL)
			{
				const DKinFitConstraint_Vertex* locVertexConstraint_CheckAgainst = dynamic_cast<const DKinFitConstraint_Vertex*>(*locIterator);
				if(locVertexConstraint_CheckAgainst == NULL)
					continue;
				if(!Check_IfKinFitResultsWillBeIdentical(static_cast<DKinFitConstraint_VertexBase*>(locVertexConstraint_ToCheck), static_cast<const DKinFitConstraint_VertexBase*>(locVertexConstraint_CheckAgainst)))
					continue;
				locMatchFoundFlag = true;
				locConstraints_CheckAgainst.erase(locIterator);
				break;
			}
			if(locSpacetimeConstraint_ToCheck != NULL)
			{
				const DKinFitConstraint_Spacetime* locSpacetimeConstraint_CheckAgainst = dynamic_cast<const DKinFitConstraint_Spacetime*>(*locIterator);
				if(locSpacetimeConstraint_CheckAgainst == NULL)
					continue;
				if(!Check_IfKinFitResultsWillBeIdentical(static_cast<DKinFitConstraint_VertexBase*>(locSpacetimeConstraint_ToCheck), static_cast<const DKinFitConstraint_VertexBase*>(locSpacetimeConstraint_CheckAgainst)))
					continue;
				if(!Check_IfKinFitResultsWillBeIdentical(locSpacetimeConstraint_ToCheck, locSpacetimeConstraint_CheckAgainst, locRFBunch_ToCheck, locRFBunch_CheckAgainst))
					continue;
				locMatchFoundFlag = true;
				locConstraints_CheckAgainst.erase(locIterator);
				break;
			}

		}
		if(!locMatchFoundFlag)
			return false;
	}
	return true;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_P4* locConstraint_ToCheck, const DKinFitConstraint_P4* locConstraint_CheckAgainst)
{
	if(!Check_IfKinFitResultsWillBeIdentical(locConstraint_ToCheck->Get_InitialParticles(), locConstraint_CheckAgainst->Get_InitialParticles()))
		return false;
	if(!Check_IfKinFitResultsWillBeIdentical(locConstraint_ToCheck->Get_FinalParticles(), locConstraint_CheckAgainst->Get_FinalParticles()))
		return false;
	return true;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_VertexBase* locConstraint_ToCheck, const DKinFitConstraint_VertexBase* locConstraint_CheckAgainst)
{
	if(!Check_IfKinFitResultsWillBeIdentical(locConstraint_ToCheck->Get_FullConstrainParticles(), locConstraint_CheckAgainst->Get_FullConstrainParticles()))
		return false;
	if(!Check_IfKinFitResultsWillBeIdentical(locConstraint_ToCheck->Get_NoConstrainParticles(), locConstraint_CheckAgainst->Get_NoConstrainParticles()))
		return false;

	deque<pair<const DKinFitParticle*, bool> > locDecayingParticles_ToCheck = locConstraint_ToCheck->Get_DecayingParticles();
	deque<pair<const DKinFitParticle*, bool> > locDecayingParticles_CheckAgainst = locConstraint_CheckAgainst->Get_DecayingParticles();
	if(locDecayingParticles_ToCheck.size() != locDecayingParticles_CheckAgainst.size())
		return false;

	//compare sources of objects //note particles can be in a different order!
	for(size_t loc_i = 0; loc_i < locDecayingParticles_ToCheck.size(); ++loc_i)
	{
		bool locMatchFoundFlag = false;
		deque<pair<const DKinFitParticle*, bool> >::iterator locIterator = locDecayingParticles_CheckAgainst.begin();
		for(; locIterator != locDecayingParticles_CheckAgainst.end(); ++locIterator)
		{
			//check pid & flag
			if((*locIterator).first->Get_PID() != locDecayingParticles_ToCheck[loc_i].first->Get_PID())
				continue;
			if((*locIterator).second != locDecayingParticles_ToCheck[loc_i].second)
				continue;

			locDecayingParticles_CheckAgainst.erase(locIterator);
			locMatchFoundFlag = true;
			break;
		}
		if(!locMatchFoundFlag)
			return false;
	}
	return true;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(DKinFitConstraint_Spacetime* locConstraint_ToCheck, const DKinFitConstraint_Spacetime* locConstraint_CheckAgainst, const DEventRFBunch* locRFBunch_ToCheck, const DEventRFBunch* locRFBunch_CheckAgainst)
{
	if(!Check_IfKinFitResultsWillBeIdentical(locConstraint_ToCheck->Get_OnlyConstrainTimeParticles(), locConstraint_CheckAgainst->Get_OnlyConstrainTimeParticles()))
		return false;
	if(locConstraint_ToCheck->Get_UseRFTimeFlag() != locConstraint_CheckAgainst->Get_UseRFTimeFlag())
		return false;
	if(locConstraint_ToCheck->Get_BeamParticle() != locConstraint_CheckAgainst->Get_BeamParticle())
		return false;
	if(locConstraint_ToCheck->Get_UseRFTimeFlag() && (locRFBunch_ToCheck != locRFBunch_CheckAgainst))
		return false;

	return true;
}

bool DKinFitResults_factory::Check_IfKinFitResultsWillBeIdentical(deque<const DKinFitParticle*> locParticles_ToCheck, deque<const DKinFitParticle*> locParticles_CheckAgainst)
{
	if(locParticles_ToCheck.size() != locParticles_CheckAgainst.size())
		return false;

	//compare sources of objects //note particles can be in a different order!
	for(size_t loc_i = 0; loc_i < locParticles_ToCheck.size(); ++loc_i)
	{
		const DKinematicData* locKinematicData_ToCheck = dKinFitter.Get_Source_FromInput(locParticles_ToCheck[loc_i]);
		bool locMatchFoundFlag = false;
		deque<const DKinFitParticle*>::iterator locIterator = locParticles_CheckAgainst.begin();
		for(; locIterator != locParticles_CheckAgainst.end(); ++locIterator)
		{
			const DKinematicData* locKinematicData_CheckAgainst = dKinFitter.Get_Source_FromInput(*locIterator);
			if(locKinematicData_ToCheck != locKinematicData_CheckAgainst)
				continue;
			if(locKinematicData_ToCheck == NULL)
			{
				//e.g. both particles missing/decaying/target: check type & pid
				if((*locIterator)->Get_PID() != locParticles_ToCheck[loc_i]->Get_PID())
					continue;
				if((*locIterator)->Get_KinFitParticleType() != locParticles_ToCheck[loc_i]->Get_KinFitParticleType())
					continue;
			}
			locParticles_CheckAgainst.erase(locIterator);
			locMatchFoundFlag = true;
			break;
		}
		if(!locMatchFoundFlag)
			return false;
	}
	return true;
}

bool DKinFitResults_factory::Create_KinFitConstraints(const DParticleCombo* locParticleCombo, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locDecayingKinFitParticles, deque<DKinFitConstraint*>& locOriginalConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints)
{
	if(dDebugLevel > 0)
		cout << "DKinFitResults_factory: Create Particles" << endl;

	locOriginalConstraints.clear();
	locSortedConstraints.clear();
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles, locFinalKinFitParticles;

	const DParticleComboStep* locParticleComboStep;
	int locDecayStepIndex;
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	const DKinematicData* locKinematicData;
	const DKinFitParticle* locKinFitParticle;
	Particle_t locPID;
	DKinFitType locKinFitType = locReaction->Get_KinFitType();

	//MAKE PARTICLES
	locInitialKinFitParticles.clear();
	locFinalKinFitParticles.clear();
	map<size_t, const DKinFitParticle*> locDecayingParticleStepMap; //the map key is the step at which the particle decays at
	bool locFirstParticleIsBeamFlag = (locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticleID() == Gamma);
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		deque<const DKinFitParticle*> locInitialKinFitParticles_Step;
		deque<const DKinFitParticle*> locFinalKinFitParticles_Step;
		locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);

		//initial particle
		locPID = locParticleComboStep->Get_InitialParticleID();
		if(locPID == Gamma)
		{
			const DBeamPhoton* locBeamPhoton = static_cast<const DBeamPhoton*>(locParticleComboStep->Get_InitialParticle());
			locKinFitParticle = dKinFitter.Make_BeamParticle(locBeamPhoton);
			locInitialKinFitParticles_Step.push_back(locKinFitParticle);
		}
		else if((!locFirstParticleIsBeamFlag) && (loc_i == 0)) //add decaying particle (production not specified, only decay)
		{
			locKinFitParticle = dKinFitter.Make_DecayingParticle(locPID);
			deque<const DKinematicData*> locFinalDecayProducts;
			locParticleCombo->Get_DecayChainParticles_Measured(loc_i, locFinalDecayProducts);
			pair<Particle_t, deque<const DKinematicData*> > locDecayInfo(locPID, locFinalDecayProducts);
			locDecayingKinFitParticles[locKinFitParticle] = locDecayInfo;
			locInitialKinFitParticles_Step.push_back(locKinFitParticle);
			locDecayingParticleStepMap[loc_i] = locKinFitParticle;
		}
		else //decaying particle, already added when it was a final particle
			locInitialKinFitParticles_Step.push_back(locDecayingParticleStepMap[loc_i]);

		//target particle
		locPID = locParticleComboStep->Get_TargetParticleID();
		if(locPID != Unknown)
		{
			locKinFitParticle = dKinFitter.Make_TargetParticle(locPID);
			locInitialKinFitParticles_Step.push_back(locKinFitParticle);
		}

		//final state particles
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
			locKinematicData = locParticleComboStep->Get_FinalParticle(loc_j);
			locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
			if(locDecayStepIndex == -1) //missing particle
			{
				locKinFitParticle = dKinFitter.Make_MissingParticle(locPID);
				locFinalKinFitParticles_Step.push_back(locKinFitParticle);
			}
			else if(locDecayStepIndex >= 0) //add decaying particle
			{
				locKinFitParticle = dKinFitter.Make_DecayingParticle(locPID);
				deque<const DKinematicData*> locFinalDecayProducts;
				locParticleCombo->Get_DecayChainParticles_Measured(locDecayStepIndex, locFinalDecayProducts);
				pair<Particle_t, deque<const DKinematicData*> > locDecayInfo(locPID, locFinalDecayProducts);
				locDecayingKinFitParticles[locKinFitParticle] = locDecayInfo;
				locFinalKinFitParticles_Step.push_back(locKinFitParticle);
				locDecayingParticleStepMap[locDecayStepIndex] = locKinFitParticle;
			}
			else if(ParticleCharge(locKinematicData->PID()) == 0) //detected neutral
			{
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
				if(locKinFitType == d_P4Fit) //make particle
					locKinFitParticle = dKinFitter.Make_DetectedParticle(locNeutralParticleHypothesis);
				else //make shower
					locKinFitParticle = dKinFitter.Make_DetectedShower(locNeutralParticleHypothesis);
				locFinalKinFitParticles_Step.push_back(locKinFitParticle);
			}
			else //detected charged track
			{
				const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locKinematicData);
				locKinFitParticle = dKinFitter.Make_DetectedParticle(locChargedTrackHypothesis);
				locFinalKinFitParticles_Step.push_back(locKinFitParticle);
			}
		}

		//add particles to deques for constraint setup
		locInitialKinFitParticles.push_back(locInitialKinFitParticles_Step);
		locFinalKinFitParticles.push_back(locFinalKinFitParticles_Step);
	}
	if(dDebugLevel > 10)
		cout << "DKinFitResults_factory: Particles Created." << endl;

	//P4: Setup constraints
	deque<DKinFitConstraint_P4*> locP4Constraints;
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Setup P4 Constraints" << endl;
		set<size_t> locStepIndicesToHandle;
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
			locStepIndicesToHandle.insert(loc_i);
		while(!locStepIndicesToHandle.empty())
		{
			deque<const DKinFitParticle*> locInitialKinFitParticles_P4, locFinalKinFitParticles_P4;
			deque<size_t> locIncludedStepIndices;

			bool locConstrainMassFlag = true; //note that the kinfitter ignores this flag when it doesn't apply (e.g. missing-mass/full-p4 constraint)
			Setup_P4Constraint(locParticleCombo, (*locStepIndicesToHandle.begin()), locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locIncludedStepIndices, locConstrainMassFlag);

			//remove steps included in the p4 constraint from the to-handle deque
			for(size_t loc_i = 0; loc_i < locIncludedStepIndices.size(); ++loc_i)
				locStepIndicesToHandle.erase(locIncludedStepIndices[loc_i]);

			//make constraint
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Create P4 Constraint" << endl;
			DKinFitConstraint_P4* locConstraint = dKinFitter.Make_P4Constraint(locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locConstrainMassFlag);
			locOriginalConstraints.push_back(dynamic_cast<DKinFitConstraint*>(locConstraint));
			locP4Constraints.push_back(locConstraint);
		}
	}

	//VERTEX OR SPACETIME: Group particles by detached vertex (one deque for each constraint/vertex)
	set<DKinFitConstraint_VertexBase*> locOriginalVertexBaseConstraints;
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Setup Vertex & Spacetime Constraints" << endl;
		set<size_t> locStepIndicesToHandle;
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
			locStepIndicesToHandle.insert(loc_i);
		while(!locStepIndicesToHandle.empty())
		{
			deque<const DKinFitParticle*> locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex;
			deque<size_t> locIncludedStepIndices;

			Setup_VertexConstraint(locParticleCombo, (*locStepIndicesToHandle.begin()), locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, locIncludedStepIndices);

			//remove steps included in the vertex constraint from the to-handle deque
			for(size_t loc_i = 0; loc_i < locIncludedStepIndices.size(); ++loc_i)
				locStepIndicesToHandle.erase(locIncludedStepIndices[loc_i]);

			if((locKinFitType == d_VertexFit) || (locKinFitType == d_P4AndVertexFit))
			{
				DKinFitConstraint_Vertex* locConstraint = dKinFitter.Make_VertexConstraint(locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, TVector3());
				//string
				string locConstraintString = string("#it{x}^{3}_{");
				for(size_t loc_j = 0; loc_j < locInitialKinFitParticles_Vertex.size(); ++loc_j)
				{
					if((loc_j == 0) || (Is_FinalStateParticle(PDGtoPType(locInitialKinFitParticles_Vertex[loc_j]->Get_PID())) == 1))
						locConstraintString += ParticleName_ROOT(PDGtoPType(locInitialKinFitParticles_Vertex[loc_j]->Get_PID()));
				}
				locConstraintString += "#rightarrow";
				for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_Vertex.size(); ++loc_j)
				{
					string locParticleString = ParticleName_ROOT(PDGtoPType(locFinalKinFitParticles_Vertex[loc_j]->Get_PID()));
					if(locFinalKinFitParticles_Vertex[loc_j]->Get_KinFitParticleType() == d_MissingParticle)
						locConstraintString += string("(") + locParticleString + string(")");
					else
						locConstraintString += locParticleString;
				}
				locConstraintString += "}";
				locConstraint->Set_ConstraintString(locConstraintString);
				//save
				locOriginalConstraints.push_back(dynamic_cast<DKinFitConstraint*>(locConstraint));
				locOriginalVertexBaseConstraints.insert(dynamic_cast<DKinFitConstraint_VertexBase*>(locConstraint));
			}
			if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
			{
				bool locUseRFTimeFlag = false;
				if(locFirstParticleIsBeamFlag && (!locInitialKinFitParticles_Vertex.empty()))
					locUseRFTimeFlag = ((locInitialKinFitParticles_Vertex[0]->Get_Charge() == 0) && (!(locInitialKinFitParticles_Vertex[0]->Get_Mass() > 0.0))); //true if photon
				DKinFitConstraint_Spacetime* locConstraint = dKinFitter.Make_SpacetimeConstraint(locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, locUseRFTimeFlag, TVector3(), 0.0);
				//string
				string locConstraintString = string("#it{x}^{4}_{");
				for(size_t loc_j = 0; loc_j < locInitialKinFitParticles_Vertex.size(); ++loc_j)
				{
					if((loc_j == 0) || (Is_FinalStateParticle(PDGtoPType(locInitialKinFitParticles_Vertex[loc_j]->Get_PID())) == 1))
						locConstraintString += ParticleName_ROOT(PDGtoPType(locInitialKinFitParticles_Vertex[loc_j]->Get_PID()));
				}
				locConstraintString += "#rightarrow";
				for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_Vertex.size(); ++loc_j)
				{
					string locParticleString = ParticleName_ROOT(PDGtoPType(locFinalKinFitParticles_Vertex[loc_j]->Get_PID()));
					if(locFinalKinFitParticles_Vertex[loc_j]->Get_KinFitParticleType() == d_MissingParticle)
						locConstraintString += string("(") + locParticleString + string(")");
					else
						locConstraintString += locParticleString;
				}
				locConstraintString += "}";
				locConstraint->Set_ConstraintString(locConstraintString);
				//save
				locOriginalConstraints.push_back(dynamic_cast<DKinFitConstraint*>(locConstraint));
				locOriginalVertexBaseConstraints.insert(dynamic_cast<DKinFitConstraint_VertexBase*>(locConstraint));
			}
		}
	}
	if(dDebugLevel > 10)
		cout << "DKinFitResults_factory: All Constraints Created." << endl;

	if(dDebugLevel > 10)
		cout << "DKinFitResults_factory: Sort Constraints." << endl;

	if(!dKinFitter.Sort_Constraints(locOriginalConstraints, locSortedConstraints))
		return false; //invalid vertex constraints are skimmed, but invalid p4 constraints cause return-false!!

	if(dDebugLevel > 10)
		cout << "DKinFitResults_factory: Constraints Sorted." << endl;

	//search for neutral showers whose vertex fits got rejected: replace them with neutral particles!
	for(size_t loc_i = 0; loc_i < locSortedConstraints.size(); ++loc_i)
		locOriginalVertexBaseConstraints.erase(locSortedConstraints[loc_i].first);
	//only the rejected constraints remain in locOriginalVertexBaseConstraints. loop over them and look for neutrals
	set<DKinFitConstraint_VertexBase*>::iterator locSetIterator = locOriginalVertexBaseConstraints.begin();
	set<const DKinFitParticle*> locNeutralShowersToReplace;
	for(; locSetIterator != locOriginalVertexBaseConstraints.end(); ++locSetIterator)
	{
		deque<const DKinFitParticle*> locNoConstrainParticles = (*locSetIterator)->Get_NoConstrainParticles();
		for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
		{
			if(locNoConstrainParticles[loc_i]->Get_IsNeutralShowerFlag()) //uh oh, need to create a particle for its p4 constraint
				locNeutralShowersToReplace.insert(locNoConstrainParticles[loc_i]);
		}
		//check spacetime constraints
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(*locSetIterator);
		if(locSpacetimeConstraint == NULL)
			continue;
		deque<const DKinFitParticle*> locOnlyConstrainTimeParticles = locSpacetimeConstraint->Get_OnlyConstrainTimeParticles();
		for(size_t loc_i = 0; loc_i < locOnlyConstrainTimeParticles.size(); ++loc_i)
		{
			if(locOnlyConstrainTimeParticles[loc_i]->Get_IsNeutralShowerFlag()) //uh oh, need to create a particle for its p4 constraint
				locNeutralShowersToReplace.insert(locNoConstrainParticles[loc_i]);
		}
	}
	if(!locNeutralShowersToReplace.empty())
	{
		//loop over p4 constraints, looking for the neutral showers we need to replace
		map<const DKinFitParticle*, const DKinematicData*> locKinematicDataMapping;
		dKinFitter.Get_ParticleMapping_InputToSource(locKinematicDataMapping);
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		{
			deque<const DKinFitParticle*> locFinalParticles = locP4Constraints[loc_i]->Get_FinalParticles();
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				if(locNeutralShowersToReplace.find(locFinalParticles[loc_j]) == locNeutralShowersToReplace.end())
					continue;

				const DKinematicData* locKinematicData = locKinematicDataMapping[locFinalParticles[loc_j]];
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
				const DKinFitParticle* locNewKinFitParticle = dKinFitter.Make_DetectedParticle(locNeutralParticleHypothesis);
				locP4Constraints[loc_i]->Replace_Particle(locFinalParticles[loc_j], false, locNewKinFitParticle); //replace!!!
			}
		}
	}

	//strip skimmed vertex/time constraints from original constraints
	deque<DKinFitConstraint*>::iterator locConstraintIterator = locOriginalConstraints.begin();
	for(; locConstraintIterator != locOriginalConstraints.end();)
	{
		if(dynamic_cast<DKinFitConstraint_VertexBase*>(*locConstraintIterator) == NULL)
		{
			++locConstraintIterator;
			continue;
		}

		//see if can find this constraint in the sorted constraints (may have been skimmed out)
		bool locConstraintFoundFlag = false;
		for(size_t loc_j = 0; loc_j < locSortedConstraints.size(); ++loc_j)
		{
			if(locSortedConstraints[loc_j].first != (*locConstraintIterator))
				continue;
			locConstraintFoundFlag = true;
			break;
		}
		if(locConstraintFoundFlag)
			++locConstraintIterator;
		else
			locConstraintIterator = locOriginalConstraints.erase(locConstraintIterator); //skimmed out, don't include
	}

	return true;
}

bool DKinFitResults_factory::Setup_KinFit(DKinFitType locKinFitType, const deque<DKinFitConstraint*>& locOriginalConstraints, const DEventRFBunch* locEventRFBunch, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints)
{
	//group p4 constraints together
	deque<DKinFitConstraint_P4*> locP4Constraints;
	for(size_t loc_i = 0; loc_i < locOriginalConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locOriginalConstraints[loc_i]);
		if(locP4Constraint != NULL)
			locP4Constraints.push_back(locP4Constraint);
	}

	//check if don't need to bother with initial vertex guesses: if so, just set constraints and return
	if((locKinFitType == d_P4Fit) || locSortedConstraints.empty())
	{
		dKinFitter.Reset_NewFit();
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
			dKinFitter.Set_Constraint(locP4Constraints[loc_i]);
		return true;
	}

	double locRFTime = locEventRFBunch->dTime;
	double locRFUncertainty = sqrt(locEventRFBunch->dTimeVariance);

	//METHOD: for each vertex (in order):
	//1) kinfit vertex only
	//2) if this is the last vertex, then exit
	//3) if possible (see 4), kinfit decay (vertex-p4) (all p4s at that vertex at once)
		//this gives decaying particles a full covariance matrix so that they can be used in subsequent vertex kinfits
	//4) if 3) is not possible because there are >= 2 locally-open-ended missing/decaying no-constrain particles (or missing particle with unknown pid): perform overall p4 fit (unless already done)
	//   (p4 only (mass constraints included where possible, all p4 constraints): get p4 cov matrices for all missing & decaying particles
	bool locP4OnlyFitPerformedFlag = false;
	map<const DKinFitParticle*, pair<const DKinFitParticle*, bool> > locReconstructedParticleMap; //map from orig particle to reconstructed particle //bool is fit good/bad (true/false)

	//last resort: decaying particles that were reconstruced in a p4-only fit that were NOT at that current vertex
		//if can't perform a vertex-p4 fit, use this info instead
	map<const DKinFitParticle*, pair<const DKinFitParticle*, bool> > locReconstructedParticleLastResortMap; //last resort map from orig particle to reconstructed particle //bool is fit good/bad (true/false)

	map<const DKinFitParticle*, const DKinFitParticle*> locReconstructedDetectedParticleMap; //map from orig decaying/missing particle to replacement "detected" particle
	for(size_t loc_i = 0; loc_i < locSortedConstraints.size(); ++loc_i)
	{
		/************************************************************** VERTEX-ONLY FIT **************************************************************/

		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Init Fit " << loc_i + 1 << " of " << locSortedConstraints.size() << ": Setup Vertex-only fit." << endl;

		DKinFitConstraint_VertexBase* locOriginalVertexBaseConstraint = locSortedConstraints[loc_i].first;
		DKinFitConstraint_Vertex* locOriginalVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(locOriginalVertexBaseConstraint);
		DKinFitConstraint_Spacetime* locOriginalSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOriginalVertexBaseConstraint);

		//may modify the original constraints for this fit
		DKinFitConstraint_VertexBase* locThisFitVertexBaseConstraint = locOriginalVertexBaseConstraint;
		DKinFitConstraint_Vertex* locThisFitVertexConstraint = locOriginalVertexConstraint;
		DKinFitConstraint_Spacetime* locThisFitSpacetimeConstraint = locOriginalSpacetimeConstraint;

		//treat previously-constrained decaying particles as detected in these constraints
		bool locConstraintClonedFlag = false;
		bool locPreviousFitFailedFlag = false;
		deque<pair<const DKinFitParticle*, bool> > locDecayingParticles = locOriginalVertexBaseConstraint->Get_DecayingParticles();
		for(size_t loc_j = 0; loc_j < locDecayingParticles.size(); ++loc_j)
		{
			map<const DKinFitParticle*, pair<const DKinFitParticle*, bool> >::iterator locMapIterator = locReconstructedParticleMap.find(locDecayingParticles[loc_j].first);
			if(locMapIterator == locReconstructedParticleMap.end())
				continue;
			if(!locMapIterator->second.second)
			{
				//previous fit for this particle failed
				if(locThisFitVertexBaseConstraint->Get_FullConstrainParticles().size() >= 2)
					continue; //but this particle isn't needed for the vertex kinematic fit: don't contaminate
				locPreviousFitFailedFlag = true; //particle is needed: don't kinfit, just return init vertex guess
				//Note: the p4 of this particle is potentially from a failed p4 fit. It should be accurate enough for a crude guess though (and better than nothing)
			}
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Reconstructed particle found, replace decaying particle." << endl;
			//replace this decaying particle with a fully-constrained "detected" particle
			if(!locConstraintClonedFlag) //don't touch the original! clone orig constraint without cloning particles
			{
				if(dDebugLevel > 15)
					cout << "DKinFitResults_factory: Clone Fit." << endl;
				if(locOriginalVertexConstraint != NULL)
				{
					locThisFitVertexConstraint = dKinFitter.Clone_KinFitConstraint_Vertex(locOriginalVertexConstraint);
					locThisFitVertexBaseConstraint = dynamic_cast<DKinFitConstraint_VertexBase*>(locThisFitVertexConstraint);
				}
				else //spacetime
				{
					locThisFitSpacetimeConstraint = dKinFitter.Clone_KinFitConstraint_Spacetime(locOriginalSpacetimeConstraint);
					locThisFitVertexBaseConstraint = dynamic_cast<DKinFitConstraint_VertexBase*>(locThisFitSpacetimeConstraint);
				}
				locConstraintClonedFlag = true;
			}
			//if not already "detected," create a new "detected" particle
			const DKinFitParticle* locReconstructedParticle = locMapIterator->second.first;
			const DKinFitParticle* locParticleToTreatAsDetected = locReconstructedParticle;
			if(locReconstructedParticle->Get_KinFitParticleType() != d_DetectedParticle)
			{
				if(dDebugLevel > 15)
					cout << "DKinFitResults_factory: Make new 'Detected' Particle." << endl;
				locParticleToTreatAsDetected = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locReconstructedParticle->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), locReconstructedParticle->Get_CovarianceMatrix());
			}
			locReconstructedDetectedParticleMap[locReconstructedParticle] = locParticleToTreatAsDetected;
			//replace decaying particle with new "detected" one
			locThisFitVertexBaseConstraint->Replace_DecayingParticle(locDecayingParticles[loc_j].first, locParticleToTreatAsDetected);
		}

		//do crude vertex guess: point on DOCA-line between the two closest (by doca) particles
		deque<const DKinFitParticle*> locFullConstrainParticles = locThisFitVertexBaseConstraint->Get_FullConstrainParticles();
		DVector3 locDVertexGuess = dAnalysisUtilities->Calc_CrudeVertex(locFullConstrainParticles);
		TVector3 locVertexGuess(locDVertexGuess.X(), locDVertexGuess.Y(), locDVertexGuess.Z());
		if(dDebugLevel > 20)
			cout << "init vertex guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;
		locThisFitVertexBaseConstraint->Set_VertexGuess(locVertexGuess);
		locOriginalVertexBaseConstraint->Set_VertexGuess(locVertexGuess);

		if((locFullConstrainParticles.size() < 2) || locPreviousFitFailedFlag)
			continue; //you failed a previous setup fit earlier, and don't have the needed information. you'll have to live with what you've got

		//setup fitter for vertex-only fit
		dKinFitter.Reset_NewFit();
		double locTimeGuess = 0.0;
		const DKinFitConstraint_VertexBase* locInitFitVertexBaseConstraint = NULL;
		if(locThisFitVertexConstraint != NULL)
			locInitFitVertexBaseConstraint = dynamic_cast<const DKinFitConstraint_VertexBase*>(dKinFitter.Set_Constraint(locThisFitVertexConstraint));
		else
		{
			locTimeGuess = Calc_TimeGuess(locThisFitSpacetimeConstraint, locDVertexGuess, locRFTime);
			locThisFitSpacetimeConstraint->Set_TimeGuess(locTimeGuess);
			locOriginalSpacetimeConstraint->Set_TimeGuess(locTimeGuess);
			locInitFitVertexBaseConstraint = dynamic_cast<const DKinFitConstraint_VertexBase*>(dKinFitter.Set_Constraint(locThisFitSpacetimeConstraint));
		}

		//perform fit
		if(dDebugLevel > 0)
			cout << "Perform init vertex guess kinematic fit" << endl;
		bool locInitVertexFitStatus = dKinFitter.Fit_Reaction();
		// it is not "necessary" for this fit to converge, because we already have a vertex/time guess (determined roughly above). 
			//thus don't return false here: if the fit fails (e.g. max #iterations), try to go on with the next iteration
		if(locInitVertexFitStatus)
		{
			//save results: constraints were cloned during fit, so need to update init guesses for the next fit
			locVertexGuess = locInitFitVertexBaseConstraint->Get_CommonVertex(); //else don't update guess!
			if(dDebugLevel > 20)
				cout << "init vertex guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;
			locThisFitVertexBaseConstraint->Set_VertexGuess(locVertexGuess);
			locOriginalVertexBaseConstraint->Set_VertexGuess(locVertexGuess);
			if(locThisFitSpacetimeConstraint != NULL)
			{
				const DKinFitConstraint_Spacetime* locInitFitSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locInitFitVertexBaseConstraint);
				locThisFitSpacetimeConstraint->Set_TimeGuess(locInitFitSpacetimeConstraint->Get_CommonTime());
				locOriginalSpacetimeConstraint->Set_TimeGuess(locInitFitSpacetimeConstraint->Get_CommonTime());
			}
		}
		map<const DKinFitParticle*, const DKinFitParticle*> locVertexFitIOMap; //input to output
		dKinFitter.Get_KinFitParticleIOMap(locVertexFitIOMap);

		/*************************************************************** CHECK P4 STATUS *************************************************************/

		if(loc_i == (locSortedConstraints.size() - 1))
			break; //have all vertex guesses, don't bother with p4

		if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit))
			continue; //no p4 fits, and decaying particles in multiple vertex/time fits with no p4 fit is illegal (will be rejected later), so don't worry about it: continue

		//see how many locally-open-ended no-constrain missing/decaying particles there are: start with decaying
			//this is an init/intermediary fit; all previously-fully-defined decaying particles should be treated as detected, and should NOT be in the decaying vector
			//note that fully-constrained decaying particles are treated as detected and thus aren't in the decaying particle vector
		deque<const DKinFitParticle*> locLocallyOpenEndedNoConstrainParticles;
		locDecayingParticles = locThisFitVertexConstraint->Get_DecayingParticles();
		for(size_t loc_j = 0; loc_j < locDecayingParticles.size(); ++loc_j)
			locLocallyOpenEndedNoConstrainParticles.push_back(locDecayingParticles[loc_i].first);

		if(locLocallyOpenEndedNoConstrainParticles.empty())
			continue; //nothing to gain by doing a vertex-p4 fit (no decaying particles to reconstruct that would be used in a future vertex constraint): just go to next fit

		//see how many locally-open-ended no-constrain missing/decaying particles there are: get missing
		deque<const DKinFitParticle*> locNoConstrainParticles = locThisFitVertexConstraint->Get_NoConstrainParticles();
		const DKinFitParticle* locVertexConstrainedMissingParticle = NULL;
		bool locUnknownMissingParticleFlag = false;
		for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
		{
			DKinFitParticleType locKinFitParticleType = locNoConstrainParticles[loc_j]->Get_KinFitParticleType();
			if(locKinFitParticleType == d_MissingParticle)
			{
				if(locNoConstrainParticles[loc_j]->Get_PID() == 0)
					locUnknownMissingParticleFlag = true;
				locVertexConstrainedMissingParticle = locVertexFitIOMap[locNoConstrainParticles[loc_j]];
				locLocallyOpenEndedNoConstrainParticles.push_back(locNoConstrainParticles[loc_j]);
			}
		}

		/*************************************************************** VERTEX-P4 FIT ***************************************************************/

		if((locLocallyOpenEndedNoConstrainParticles.size() < 2) && (!locUnknownMissingParticleFlag))
		{
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Init Fit " << loc_i + 1 << " of " << locSortedConstraints.size() << ": Setup Vertex-P4 fit." << endl;

			//perform a vertex-p4 fit to get: better vertex guess & full cov matrix of remaining missing/decaying particle
			dKinFitter.Reset_NewFit();
			const DKinFitConstraint_VertexBase* locP4VertexFitVertexBaseConstraint = NULL;
			if(locThisFitVertexConstraint != NULL)
				locP4VertexFitVertexBaseConstraint = dynamic_cast<const DKinFitConstraint_VertexBase*>(dKinFitter.Set_Constraint(locThisFitVertexConstraint));
			else
				locP4VertexFitVertexBaseConstraint = dynamic_cast<const DKinFitConstraint_VertexBase*>(dKinFitter.Set_Constraint(locThisFitSpacetimeConstraint));

			//get & potentially modify p4 constraints: treat previously-constrained decaying/missing particles as detected in these constraints
			set<DKinFitConstraint_P4*> locTempP4Constraints = locSortedConstraints[loc_i].second;
			for(set<DKinFitConstraint_P4*>::iterator locP4Iterator = locTempP4Constraints.begin(); locP4Iterator != locTempP4Constraints.end(); ++locP4Iterator)
			{
				DKinFitConstraint_P4* locOriginalP4Constraint = *locP4Iterator;
				//get missing/decaying particles in this constraint
				deque<pair<const DKinFitParticle*, bool> > locMissingOrDecayingParticles; //bool is true/false if in initial/final state
				deque<const DKinFitParticle*> locInitialParticles = locOriginalP4Constraint->Get_InitialParticles();
				for(size_t loc_j = 0; loc_j < locInitialParticles.size(); ++loc_j)
				{
					DKinFitParticleType locKinFitParticleType = locInitialParticles[loc_j]->Get_KinFitParticleType();
					if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle))
						locMissingOrDecayingParticles.push_back(pair<const DKinFitParticle*, bool>(locInitialParticles[loc_j], true));
				}
				deque<const DKinFitParticle*> locFinalParticles = locOriginalP4Constraint->Get_FinalParticles();
				for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
				{
					DKinFitParticleType locKinFitParticleType = locFinalParticles[loc_j]->Get_KinFitParticleType();
					if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle))
						locMissingOrDecayingParticles.push_back(pair<const DKinFitParticle*, bool>(locFinalParticles[loc_j], false));
				}

				//loop over missing/decaying particles, replacing with reconstructed particles as available
				DKinFitConstraint_P4* locClonedP4Constraint = NULL;
				for(size_t loc_j = 0; loc_j < locMissingOrDecayingParticles.size(); ++loc_j)
				{
					map<const DKinFitParticle*, pair<const DKinFitParticle*, bool> >::iterator locMapIterator = locReconstructedParticleMap.find(locMissingOrDecayingParticles[loc_j].first);
					if(locMapIterator == locReconstructedParticleMap.end())
						continue;
					if(!locMapIterator->second.second)
					{
						locPreviousFitFailedFlag = true;
						break; //this particle failed in a previous fit: cannot trust results. do not perform vertex-p4 fit
					}
					//replace this decaying/missing particle with a fully-constrained "detected" particle

					//don't touch the original constraint! clone orig constraint (if haven't already) without cloning particles
					if(locClonedP4Constraint == NULL)
						locClonedP4Constraint = dKinFitter.Clone_KinFitConstraint_P4(locOriginalP4Constraint);

					const DKinFitParticle* locParticleToTreatAsDetected = NULL;
					map<const DKinFitParticle*, const DKinFitParticle*>::iterator locDetectedIterator = locReconstructedDetectedParticleMap.find(locMissingOrDecayingParticles[loc_j].first);
					if(locDetectedIterator != locReconstructedDetectedParticleMap.end())
						locParticleToTreatAsDetected = locDetectedIterator->second; //"detected" particle already made: use it
					else
					{
						//create new "detected" particle
						const DKinFitParticle* locReconstructedParticle = locMapIterator->second.first;
						if(locInitVertexFitStatus)
						{
							if(locReconstructedParticle->Get_KinFitParticleType() == d_MissingParticle)
							{
								//combine info from above vertex fit & previous p3 fit of missing particle
								//Note that if it is charged and in a B-field, then the missing momentum will be propagated between this input vertex position and the soon-to-be-kinfit one
									//However, this distance is likely to be small, and the missing p4 will likely not change over this distance
								TMatrixDSym locCovarianceMatrix = *(locVertexConstrainedMissingParticle->Get_CovarianceMatrix()); //contains v3 & t from current vertex fit
								const TMatrixDSym* locP3CovarianceMatrix = locReconstructedParticle->Get_CovarianceMatrix(); //contains p3 from overall-p4 fit
								//update with p3 from overall-p4 fit
								for(unsigned int loc_q = 0; loc_q < 3; ++loc_q)
								{
									for(unsigned int loc_r = 0; loc_r < 3; ++loc_r)
										locCovarianceMatrix(loc_q, loc_r) = (*locP3CovarianceMatrix)(loc_q, loc_r);
								}
								locParticleToTreatAsDetected = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locVertexConstrainedMissingParticle->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
							}
							else if(locReconstructedParticle->Get_KinFitParticleType() == d_DecayingParticle)
								locParticleToTreatAsDetected = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locReconstructedParticle->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), locReconstructedParticle->Get_CovarianceMatrix());
							else //already reconstructed as detected, just copy it
								locParticleToTreatAsDetected = locReconstructedParticle;
						}
						else //cannot trust init vertex/time fit
						{
							TMatrixDSym locCovarianceMatrix = *(locReconstructedParticle->Get_CovarianceMatrix()); //contains p3 from overall-p4 fit
							//force ~infinite uncertainty on vertex params //it's only a setup fit anyway
							for(unsigned int loc_q = 0; loc_q < 4; ++loc_q)
							{
								for(unsigned int loc_r = 0; loc_r < 4; ++loc_r)
									locCovarianceMatrix(loc_q + 3, loc_r + 3) = 9.9E99;
							}

							if(locReconstructedParticle->Get_KinFitParticleType() == d_MissingParticle)
								locParticleToTreatAsDetected = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), TLorentzVector(locVertexGuess, locTimeGuess), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
							else if(locReconstructedParticle->Get_KinFitParticleType() == d_DecayingParticle)
								locParticleToTreatAsDetected = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), TLorentzVector(locVertexGuess, locTimeGuess), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
							else //already reconstructed as detected, just copy it
								locParticleToTreatAsDetected = locReconstructedParticle;
						}
						locReconstructedDetectedParticleMap[locReconstructedParticle] = locParticleToTreatAsDetected;
					}

					//replace missing/decaying particle with new "detected" one
					locClonedP4Constraint->Replace_Particle(locMissingOrDecayingParticles[loc_j].first, locMissingOrDecayingParticles[loc_j].second, locParticleToTreatAsDetected);
				}
				if(locPreviousFitFailedFlag)
					break;
				if(locClonedP4Constraint != NULL)
					dKinFitter.Set_Constraint(locClonedP4Constraint);
				else
					dKinFitter.Set_Constraint(locOriginalP4Constraint);
			}

			//perform fit
			if(dDebugLevel > 0)
				cout << "Perform init p4-vertex guess kinematic fit" << endl;

			bool locP4VertexFitStatus = (!locPreviousFitFailedFlag) ? dKinFitter.Fit_Reaction() : false;
			if(locP4VertexFitStatus)
			{
				TVector3 locVertexGuess = locP4VertexFitVertexBaseConstraint->Get_CommonVertex();
				if(dDebugLevel > 20)
					cout << "init vertex guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;

				//save v3/t results
				locOriginalVertexBaseConstraint->Set_VertexGuess(locVertexGuess);
				if(locOriginalSpacetimeConstraint != NULL)
				{
					const DKinFitConstraint_Spacetime* locInitFitSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locP4VertexFitVertexBaseConstraint);
					locOriginalSpacetimeConstraint->Set_TimeGuess(locInitFitSpacetimeConstraint->Get_CommonTime());
				}

				//save reconstructed decaying/missing particles
					//assumes all particles reconstructed during vertex-p4 fit were in both p4 & vertex (as no-constrain particles) fits
				deque<const DKinFitParticle*> locNoConstrainParticles = locP4VertexFitVertexBaseConstraint->Get_NoConstrainParticles();
				for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
				{
					DKinFitParticleType locKinFitParticleType = locNoConstrainParticles[loc_j]->Get_KinFitParticleType();
					if(locKinFitParticleType != d_DecayingParticle)
						continue;
					const DKinFitParticle* locOriginalKinFitParticle = dKinFitter.Get_InputKinFitParticle(locNoConstrainParticles[loc_j]);
					locReconstructedParticleMap[locOriginalKinFitParticle] = pair<const DKinFitParticle*, bool>(locNoConstrainParticles[loc_j], true);
					break;
				}
				continue;
			}
			//if vertex-p4 fit did not converge, move onto full-p4 fit to get decaying particle info
		}

		/************************************************************* FULL, P4-ONLY FIT *************************************************************/

		//too many unknown particles: cannot do local vertex-p4 fit    OR    vertex-p4 fit failed
		if(!locP4OnlyFitPerformedFlag)
		{
			//if not done already, perform a full, p4-only fit (mass constraints are applied as possible: especially important if inclusive fit)
			//save the resulting, reconstructed decaying/missing particles for use in later fits
				//will not use decaying particles reconstructed at other vertices, except as a last resort

			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Init Fit " << loc_i + 1 << " of " << locSortedConstraints.size() << ": Setup P4-Only fit." << endl;

			//setup fit
			dKinFitter.Reset_NewFit();
			deque<const DKinFitConstraint_P4*> locInitFitP4Constraints;

			//must treat any neutral showers as neutral particles (no vertex fit here)
			map<const DKinFitParticle*, const DKinematicData*> locKinematicDataMapping;
			dKinFitter.Get_ParticleMapping_InputToSource(locKinematicDataMapping);
			for(size_t loc_j = 0; loc_j < locP4Constraints.size(); ++loc_j)
			{
				bool locConstraintClonedFlag = false;
				DKinFitConstraint_P4* locThisP4Constraint = locP4Constraints[loc_j];
				for(size_t loc_k = 0; loc_k < locP4Constraints[loc_j]->Get_FinalParticles().size(); ++loc_k)
				{
					const DKinFitParticle* locKinFitParticle = (locP4Constraints[loc_j]->Get_FinalParticles())[loc_k];
					if(!locKinFitParticle->Get_IsNeutralShowerFlag())
						continue;
					const DKinematicData* locKinematicData = locKinematicDataMapping[locKinFitParticle];
					const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
					const DKinFitParticle* locNewKinFitParticle = dKinFitter.Make_DetectedParticle(locNeutralParticleHypothesis);
					if(!locConstraintClonedFlag) //don't touch the original! clone orig constraint without cloning particles
					{
						if(dDebugLevel > 15)
							cout << "DKinFitResults_factory: Clone P4 Fit." << endl;
						locThisP4Constraint = dKinFitter.Clone_KinFitConstraint_P4(locThisP4Constraint);
					}
					locThisP4Constraint->Replace_Particle(locKinFitParticle, false, locNewKinFitParticle);
				}
				locInitFitP4Constraints.push_back(dKinFitter.Set_Constraint(locThisP4Constraint));
			}

			//do fit
			if(dDebugLevel > 0)
				cout << "Perform init full-p4 kinematic fit" << endl;
			bool locP4OnlyFitStatus = dKinFitter.Fit_Reaction();

			//save results: save all reconstructed decaying/missing particle info (if not saved already): p3 & cov p3
				//save the missing & decaying-no-constrain-particles-at-this-vertex in the reconstruction map
				//save decaying particles at other vertices TO A SPECIAL LAST RESORT MAP

			//loop over created p4 constraints, grabbing the constrained particle (decaying/missing) in each one (and saving the results)
			deque<const DKinFitParticle*> locNoVertexConstrainParticles = locInitFitVertexBaseConstraint->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locInitFitP4Constraints.size(); ++loc_j)
			{
				//even if the fit did not converge, locReconstructedParticle should have a semi-good guess for the p4
				const DKinFitParticle* locReconstructedParticle = locInitFitP4Constraints[loc_j]->Get_ConstrainedP4Particle();
				const DKinFitParticle* locOriginalParticle = dKinFitter.Get_InputKinFitParticle(locReconstructedParticle);
				if((locReconstructedParticle == NULL) || (locOriginalParticle == NULL))
					continue; //e.g. no missing particle
				//see if a reconstructed version of this particle has already been saved
				if(locReconstructedParticleMap.find(locOriginalParticle) != locReconstructedParticleMap.end())
					continue; //have already saved a reconstructed version of this particle, which is better than this new one. skip it

				//if decaying particle defined at this vertex: save the v3 & v3 cov matrix
				const DKinFitParticle* locVertexFitResult = NULL;
				for(size_t loc_k = 0; loc_k < locNoVertexConstrainParticles.size(); ++loc_k)
				{
					//need the original
					const DKinFitParticle* locVertexFitOriginal = NULL;
					map<const DKinFitParticle*, const DKinFitParticle*>::const_iterator locIOMapIterator = locVertexFitIOMap.begin();
					for(; locIOMapIterator != locVertexFitIOMap.end(); ++locIOMapIterator)
					{
						if(locIOMapIterator->second == locNoVertexConstrainParticles[loc_k])
							locVertexFitOriginal = locIOMapIterator->first;
					}
					if(locVertexFitOriginal == NULL)
						continue; //shouldn't be possible ...
					//compare originals
					if(locVertexFitOriginal != locOriginalParticle)
						continue;
					locVertexFitResult = locNoVertexConstrainParticles[loc_k];
					break;
				}

				//save the result: missing particle
				if(locOriginalParticle->Get_KinFitParticleType() != d_DecayingParticle)
				{
					if(!locP4OnlyFitStatus)
					{
						//need a non-null covariance matrix
						TMatrixDSym locCovarianceMatrix(7, 7);
						locCovarianceMatrix.ResizeTo(7, 7);
						locReconstructedParticle = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locReconstructedParticle->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
					}
					locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locReconstructedParticle, locP4OnlyFitStatus);
					continue;
				}
				else if(locVertexFitResult == NULL) //decaying particle but not in this vertex fit
				{
					if(!locP4OnlyFitStatus)
					{
						//need a non-null covariance matrix
						TMatrixDSym locCovarianceMatrix(7, 7);
						locCovarianceMatrix.ResizeTo(7, 7);
						locReconstructedParticle = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locReconstructedParticle->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
					}
					if(dDebugLevel > 10)
						cout << "DKinFitResults_factory: Save last resort particle for PID = " << locReconstructedParticle->Get_PID() << endl;
					locReconstructedParticleLastResortMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locReconstructedParticle, locP4OnlyFitStatus);
				}
				else //decaying particle in this vertex fit
				{
					//combine info from above vertex fit with this p3 fit result: create new "detected" particle
					if(locInitVertexFitStatus)
					{
						//fit converged
						TMatrixDSym locCovarianceMatrix = *(locVertexFitResult->Get_CovarianceMatrix()); //contains v3 & t from current vertex fit
						if(locP4OnlyFitStatus)
						{
							const TMatrixDSym* locP3CovarianceMatrix = locReconstructedParticle->Get_CovarianceMatrix(); //contains p3 from overall-p4 fit
							//update with p3 from overall-p4 fit
							for(unsigned int loc_q = 0; loc_q < 3; ++loc_q)
							{
								for(unsigned int loc_r = 0; loc_r < 3; ++loc_r)
									locCovarianceMatrix(loc_q, loc_r) = (*locP3CovarianceMatrix)(loc_q, loc_r);
							}
						}
						const DKinFitParticle* locNewReconstructedParticle = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), locVertexFitResult->Get_SpacetimeVertex(), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
						locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locNewReconstructedParticle, locP4OnlyFitStatus);
					}
					else
					{
						//cannot trust init vertex/time fit
						TMatrixDSym locCovarianceMatrix(7, 7);
						locCovarianceMatrix.ResizeTo(7, 7);
						//force ~infinite uncertainty on vertex params //it's only a setup fit anyway
						const DKinFitParticle* locNewReconstructedParticle = dKinFitter.Make_DetectedParticle(locReconstructedParticle->Get_PID(), locReconstructedParticle->Get_Charge(), locReconstructedParticle->Get_Mass(), TLorentzVector(locVertexGuess, locTimeGuess), locReconstructedParticle->Get_Momentum(), &locCovarianceMatrix);
						locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locNewReconstructedParticle, false); //fit failed
					}
				}
			}
			locP4OnlyFitPerformedFlag = true;
			continue;
		}

		/**************************************************** USE LAST-RESORT P4-ONLY-FIT RESULTS ****************************************************/

		//p4-only fit already performed, so instead, grab info for previously-fit decaying particles from last-resort reconstruction map
			//they weren't used earlier, because we were hoping to do a vertex-p4 fit to get those results
		//get needed decaying particles from last resort map, update their info with the latest vertex, push to the reconstruction map, and continue (don't do any fits)

		for(size_t loc_j = 0; loc_i < locLocallyOpenEndedNoConstrainParticles.size(); ++loc_j)
		{
			//get from last-resort map
			const DKinFitParticle* locOriginalParticle = locLocallyOpenEndedNoConstrainParticles[loc_j];
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Get last resort particle for PID = " << locOriginalParticle->Get_PID() << endl;

			map<const DKinFitParticle*, pair<const DKinFitParticle*, bool> >::iterator locLastResortIterator = locReconstructedParticleLastResortMap.find(locOriginalParticle);
			if(locLastResortIterator == locReconstructedParticleLastResortMap.end())
				continue; //e.g. missing

			//combine info from above vertex fit with this p3 fit result: create new "detected" particle
			const DKinFitParticle* locLastResortParticle = locLastResortIterator->second.first;
			if(!locInitVertexFitStatus)
			{
				//can't improve, use as is
				locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locLastResortParticle, false);
				continue;
			}

			const DKinFitParticle* locVertexFitResult = locVertexFitIOMap[locOriginalParticle];
			if(!locLastResortIterator->second.second)
			{
				//p4 fit failed, but can at least set vertex position
				const DKinFitParticle* locNewReconstructedParticle = dKinFitter.Make_DetectedParticle(locLastResortParticle->Get_PID(), locLastResortParticle->Get_Charge(), locLastResortParticle->Get_Mass(), locVertexFitResult->Get_SpacetimeVertex(), locLastResortParticle->Get_Momentum(), locVertexFitResult->Get_CovarianceMatrix());
				locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locNewReconstructedParticle, false);
				continue;
			}
			//p4 & vertex fit successful: combine info into new particle
			TMatrixDSym locCovarianceMatrix = *(locVertexFitResult->Get_CovarianceMatrix()); //contains v3 & t from current vertex fit
			const TMatrixDSym* locP3CovarianceMatrix = locLastResortParticle->Get_CovarianceMatrix(); //contains p3 from overall-p4 fit
			//update with p3 from overall-p4 fit
			for(unsigned int loc_q = 0; loc_q < 3; ++loc_q)
			{
				for(unsigned int loc_r = 0; loc_r < 3; ++loc_r)
					locCovarianceMatrix(loc_q, loc_r) = (*locP3CovarianceMatrix)(loc_q, loc_r);
			}
			const DKinFitParticle* locNewReconstructedParticle = dKinFitter.Make_DetectedParticle(locLastResortParticle->Get_PID(), locLastResortParticle->Get_Charge(), locLastResortParticle->Get_Mass(), locVertexFitResult->Get_SpacetimeVertex(), locLastResortParticle->Get_Momentum(), &locCovarianceMatrix);
			locReconstructedParticleMap[locOriginalParticle] = pair<const DKinFitParticle*, bool>(locNewReconstructedParticle, true);
		}
	}

	if(dDebugLevel > 10)
		cout << "DKinFitResults_factory: Init Guesses found; set all constraints for master fit." << endl;

	dKinFitter.Reset_NewFit();

	//ADD CONSTRAINTS
	for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		dKinFitter.Set_Constraint(locP4Constraints[loc_i]);

	//use locOriginalConstraints to keep constraints in same order (doesn't really matter, but confidence level histogram looks better this way)
	const DKinFitParticle* locBeamParticleForRF = NULL;
	for(size_t loc_i = 0; loc_i < locOriginalConstraints.size(); ++loc_i)
	{
		if(dynamic_cast<DKinFitConstraint_VertexBase*>(locOriginalConstraints[loc_i]) == NULL)
			continue;

		DKinFitConstraint_Vertex* locVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(locOriginalConstraints[loc_i]);
		if(locVertexConstraint != NULL)
		{
			dKinFitter.Set_Constraint(locVertexConstraint);
			if(dDebugLevel > 10)
			{
				TVector3 locVertex = locVertexConstraint->Get_CommonVertex();
				cout << "final vertex guess (xyz) = " << locVertex.X() << ", " << locVertex.Y() << ", " << locVertex.Z() << endl;
			}
			continue;
		}
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOriginalConstraints[loc_i]);
		if(locSpacetimeConstraint != NULL)
		{
			if((locBeamParticleForRF == NULL) && locSpacetimeConstraint->Get_UseRFTimeFlag())
				locBeamParticleForRF = locSpacetimeConstraint->Get_BeamParticle();
			dKinFitter.Set_Constraint(locSpacetimeConstraint);
			if(dDebugLevel > 10)
			{
				TLorentzVector locLorentzVector = locSpacetimeConstraint->Get_CommonSpacetimeVertex();
				cout << "final spacetime guess (xyzt) = " << locLorentzVector.X() << ", " << locLorentzVector.Y() << ", " << locLorentzVector.Z() << ", " << locLorentzVector.T() << endl;
			}
		}
	}

	//Add RF time if needed
	if(locBeamParticleForRF != NULL)
		dKinFitter.Set_RFTime(locRFTime, locRFUncertainty, locBeamParticleForRF);

	return true;
}

void DKinFitResults_factory::Setup_P4Constraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_P4, deque<const DKinFitParticle*>& locFinalKinFitParticles_P4, deque<size_t>& locIncludedStepIndices, bool& locConstrainMassFlag)
{
	locIncludedStepIndices.push_back(locStepIndex);
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	Particle_t locPID;

	//initial particle
	locPID = locParticleComboStep->Get_InitialParticleID();
	if(locPID == Gamma)
		locInitialKinFitParticles_P4.push_back(locInitialKinFitParticles[locStepIndex][0]);
	else //decaying particle
	{
		locConstrainMassFlag = locParticleCombo->Get_ApplyKinFitMassConstraintOnInitialParticleFlag(locStepIndex);
		if(IsFixedMass(locPID)) //else is decaying resonance particle: don't constrain (don't recurse either (already done!))
			locInitialKinFitParticles_P4.push_back(locInitialKinFitParticles[locStepIndex][0]);
	}

	//target particle
	if(locParticleComboStep->Get_TargetParticleID() != Unknown)
		locInitialKinFitParticles_P4.push_back(locInitialKinFitParticles[locStepIndex][1]);

	//final state particles
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
		if(locDecayStepIndex == -1) //missing particle
			locFinalKinFitParticles_P4.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
		else if(locDecayStepIndex >= 0) //decaying particle
		{
			if(!IsFixedMass(locPID)) //go to the next step!!
				Setup_P4Constraint(locParticleCombo, locDecayStepIndex, locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locIncludedStepIndices, locConstrainMassFlag);
			else
				locFinalKinFitParticles_P4.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
		}
		else //detected particle or shower
			locFinalKinFitParticles_P4.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
	}
}

void DKinFitResults_factory::Setup_VertexConstraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_Vertex, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex, deque<size_t>& locIncludedStepIndices)
{
	locIncludedStepIndices.push_back(locStepIndex);
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	Particle_t locPID;

	//initial particle
	locPID = locParticleComboStep->Get_InitialParticleID();
	if(locPID == Gamma)
		locInitialKinFitParticles_Vertex.push_back(locInitialKinFitParticles[locStepIndex][0]);
	else //decaying particle
		locInitialKinFitParticles_Vertex.push_back(locInitialKinFitParticles[locStepIndex][0]);

	//target particle
	if(locParticleComboStep->Get_TargetParticleID() != Unknown)
		locInitialKinFitParticles_Vertex.push_back(locInitialKinFitParticles[locStepIndex][1]);

	//final state particles
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
		if(locDecayStepIndex == -1) //missing particle
			locFinalKinFitParticles_Vertex.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
		else if(locDecayStepIndex >= 0) //decaying particle
		{
			if(IsDetachedVertex(locPID))
			{
				if(dLinkVerticesFlag)
					locFinalKinFitParticles_Vertex.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
			}
			else //go to the next step!!
				Setup_VertexConstraint(locParticleCombo, locDecayStepIndex, locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, locIncludedStepIndices);
		}
		else //detected particle or shower
			locFinalKinFitParticles_Vertex.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
	}
}

double DKinFitResults_factory::Calc_TimeGuess(const DKinFitConstraint_Spacetime* locConstraint, DVector3 locVertexGuess, double locRFTime)
{
	//if RF: propagate RF time to the vertex-z of the init vertex guess
	bool locUseRFTimeFlag = locConstraint->Get_UseRFTimeFlag();
	if(locUseRFTimeFlag)
		return (locRFTime + (locVertexGuess.Z() - dTargetZCenter)/29.9792458);

	//propagate each track time to the DOCA to the init vertex guess and average them
	deque<const DKinFitParticle*> locTimeFindParticles = locConstraint->Get_FullConstrainParticles();
	deque<const DKinFitParticle*> locOnlyTimeFindParticles = locConstraint->Get_OnlyConstrainTimeParticles();
	locTimeFindParticles.insert(locTimeFindParticles.end(), locOnlyTimeFindParticles.begin(), locOnlyTimeFindParticles.end());
	return dAnalysisUtilities->Calc_CrudeTime(locTimeFindParticles, locVertexGuess);
}

void DKinFitResults_factory::Build_KinFitResults(const DParticleCombo* locParticleCombo, const map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locInitDecayingKinFitParticles, deque<DKinFitConstraint*>& locOriginalConstraints)
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	DKinFitResults* locKinFitResults = new DKinFitResults();
	locKinFitResults->Add_ParticleCombo(locParticleCombo);

	locKinFitResults->Set_KinFitType(locReaction->Get_KinFitType());
	locKinFitResults->Set_ConfidenceLevel(dKinFitter.Get_ConfidenceLevel());
	locKinFitResults->Set_ChiSq(dKinFitter.Get_ChiSq());
	locKinFitResults->Set_NDF(dKinFitter.Get_NDF());
	locKinFitResults->Set_VEta(dKinFitter.Get_VEta());
	locKinFitResults->Set_VXi(dKinFitter.Get_VXi());
	locKinFitResults->Set_V(dKinFitter.Get_V());
	locKinFitResults->Set_NumConstraints(dKinFitter.Get_NumConstraintEquations());
	locKinFitResults->Set_NumUnknowns(dKinFitter.Get_NumUnknowns());
	locKinFitResults->Set_OriginalKinFitConstraints(locOriginalConstraints);
	locKinFitResults->Set_InputDecayingParticleInfo(locInitDecayingKinFitParticles);

	//save the constraints
	deque<const DKinFitConstraint*> locKinFitConstraints;
	dKinFitter.Get_KinFitConstraints(locKinFitConstraints);

	string locBeamTargetString = "";
	if(locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticleID() == Gamma)
		locBeamTargetString += ParticleName_ROOT(Gamma);
	if(locParticleCombo->Get_ParticleComboStep(0)->Get_TargetParticleID() != Unknown)
		locBeamTargetString += ParticleName_ROOT(locParticleCombo->Get_ParticleComboStep(0)->Get_TargetParticleID());
	locBeamTargetString += "#rightarrow";

	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		const DKinFitConstraint_P4* locP4Constraint = dynamic_cast<const DKinFitConstraint_P4*>(locKinFitConstraints[loc_i]);
		if(locP4Constraint != NULL)
		{
			string locConstraintString = "";
			if(locP4Constraint->Get_IsActualP4ConstraintFlag())
				locConstraintString = "#it{p}^{4}";
			else if(locP4Constraint->Get_ConstrainInitialParticleMassFlag())
			{
				const DKinFitParticle* locKinFitParticle = (locP4Constraint->Get_InitialParticles())[0];
				locConstraintString = string("#it{m}_{") + string(ParticleName_ROOT(PDGtoPType(locKinFitParticle->Get_PID()))) + string("}");
			}
			(const_cast<DKinFitConstraint_P4*>(locP4Constraint))->Set_ConstraintString(locConstraintString);
		}
		else //vertex/time constraint, set colors based on if used to constrain
		{
			const DKinFitConstraint_VertexBase* locVertexBaseConstraint = dynamic_cast<const DKinFitConstraint_VertexBase*>(locKinFitConstraints[loc_i]);
			TString locConstraintTString = locVertexBaseConstraint->Get_ConstraintString();
			deque<const DKinFitParticle*> locFullConstrainParticles = locVertexBaseConstraint->Get_FullConstrainParticles();
			set<int> locPIDsReplaced;
			for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
			{
				int locPID = locFullConstrainParticles[loc_j]->Get_PID();
				if(locFullConstrainParticles[loc_j]->Get_KinFitParticleType() == d_BeamParticle)
				{
					string locParticleString = ParticleName_ROOT(PDGtoPType(locPID));
					locConstraintTString.Replace(0, locParticleString.size(), (string("#color[4]{") + locParticleString + string("}")).c_str());
					locBeamTargetString.insert(locParticleString.size(), "}");
					continue;
				}
				if(locPIDsReplaced.find(locPID) != locPIDsReplaced.end())
					continue;
				string locParticleString = ParticleName_ROOT(PDGtoPType(locPID));
				string locMissingParticleString = string("(") + locParticleString + string(")");
				locConstraintTString.ReplaceAll(locMissingParticleString.c_str(), "__TEMPORARY__"); //temporary //don't want to make missing particles blue!
				locConstraintTString.ReplaceAll(locBeamTargetString.c_str(), "__BEAMTARGET__"); //temporary //don't want to make beam/target particles blue! (here)

				string locNewParticleString = string("#color[4]{") + locParticleString + string("}"); //blue
				locConstraintTString.ReplaceAll(locParticleString.c_str(), locNewParticleString.c_str());
				locConstraintTString.ReplaceAll("__TEMPORARY__", locMissingParticleString.c_str()); //restore
				locConstraintTString.ReplaceAll("__BEAMTARGET__", locBeamTargetString.c_str()); //restore
				locPIDsReplaced.insert(locPID);
			}
			const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locVertexBaseConstraint);
			if(locSpacetimeConstraint != NULL)
			{
				deque<const DKinFitParticle*> locOnlyConstrainTimeParticles = locSpacetimeConstraint->Get_OnlyConstrainTimeParticles();
				for(size_t loc_j = 0; loc_j < locOnlyConstrainTimeParticles.size(); ++loc_j)
				{
					int locPID = locOnlyConstrainTimeParticles[loc_j]->Get_PID();
					if(locOnlyConstrainTimeParticles[loc_j]->Get_KinFitParticleType() == d_BeamParticle)
					{
						string locParticleString = ParticleName_ROOT(PDGtoPType(locPID));
						locConstraintTString.Replace(0, locParticleString.size(), (string("#color[4]{") + locParticleString + string("}")).c_str());
						locBeamTargetString.insert(locParticleString.size(), "}");
						continue;
					}
					if(locPIDsReplaced.find(locPID) != locPIDsReplaced.end())
						continue;
					string locParticleString = ParticleName_ROOT(PDGtoPType(locPID));
					string locMissingParticleString = string("(") + locParticleString + string(")");
					locConstraintTString.ReplaceAll(locMissingParticleString.c_str(), "__TEMPORARY__"); //temporary //don't want to make missing particles blue!
					locConstraintTString.ReplaceAll(locBeamTargetString.c_str(), "__BEAMTARGET__"); //temporary //don't want to make beam/target particles blue! (here)

					string locNewParticleString = string("#color[4]{") + locParticleString + string("}"); //blue
					locConstraintTString.ReplaceAll(locParticleString.c_str(), locNewParticleString.c_str());
					locConstraintTString.ReplaceAll("__TEMPORARY__", locMissingParticleString.c_str()); //restore
					locConstraintTString.ReplaceAll("__BEAMTARGET__", locBeamTargetString.c_str()); //restore
					locPIDsReplaced.insert(locPID);
				}
			}
			(const_cast<DKinFitConstraint_VertexBase*>(locVertexBaseConstraint))->Set_ConstraintString((const char*)locConstraintTString);
		}
	}

	locKinFitResults->Set_KinFitConstraints(locKinFitConstraints);

	//map of particle data
	map<const DKinFitParticle*, const DKinematicData*> locParticleMapping_Output;
	dKinFitter.Get_ParticleMapping_OutputToSource(locParticleMapping_Output);
	locKinFitResults->Set_ParticleMapping(locParticleMapping_Output);

	//reverse map of particle data (& missing)
	map<const DKinematicData*, const DKinFitParticle*> locReverseParticleMapping;
	for(map<const DKinFitParticle*, const DKinematicData*>::iterator locIterator = locParticleMapping_Output.begin(); locIterator != locParticleMapping_Output.end(); ++locIterator)
	{
		if(locIterator->first->Get_KinFitParticleType() == d_MissingParticle)
			locKinFitResults->Set_MissingParticle(locIterator->first); //missing
		if(locIterator->second == NULL)
			continue; //not detected (e.g. missing, decaying, target)
		locReverseParticleMapping[locIterator->second] = locIterator->first;
	}
	locKinFitResults->Set_ReverseParticleMapping(locReverseParticleMapping);

	//decaying particles //input to function is the constructed decaying particles; must save the final decaying particles
	map<const DKinFitParticle*, const DKinFitParticle*> locKinFitParticleIOMap; //map from constructed-kinfit-particle to final-kinfit-particle
	dKinFitter.Get_KinFitParticleIOMap(locKinFitParticleIOMap);
	for(map<const DKinFitParticle*, const DKinFitParticle*>::iterator locIterator = locKinFitParticleIOMap.begin(); locIterator != locKinFitParticleIOMap.end(); ++locIterator)
	{
		map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >::const_iterator locDecayingIterator = locInitDecayingKinFitParticles.find(locIterator->first);
		if(locDecayingIterator == locInitDecayingKinFitParticles.end())
			continue; //not a decaying particle
		locKinFitResults->Add_DecayingParticle(locDecayingIterator->second, locKinFitParticleIOMap[locDecayingIterator->first]);
	}

	//pulls
	map<const DKinematicData*, map<DKinFitPullType, double> > locPulls;
	dKinFitter.Get_Pulls(locPulls);
	locKinFitResults->Set_Pulls(locPulls);

	_data.push_back(locKinFitResults);
}

//------------------
// erun
//------------------
jerror_t DKinFitResults_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DKinFitResults_factory::fini(void)
{
	return NOERROR;
}

