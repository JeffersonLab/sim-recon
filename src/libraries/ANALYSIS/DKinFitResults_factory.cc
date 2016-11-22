#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DKinFitResults_factory.h"

//------------------
// init
//------------------
jerror_t DKinFitResults_factory::init(void)
{
	dKinFitDebugLevel = 0;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DKinFitResults_factory::brun(jana::JEventLoop* locEventLoop, int32_t runnumber)
{
	//CREATE FIT UTILS AND FITTER
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);
	dKinFitter = new DKinFitter(dKinFitUtils);

	gPARMS->SetDefaultParameter("KINFIT:DEBUGLEVEL", dKinFitDebugLevel);
	dKinFitter->Set_DebugLevel(dKinFitDebugLevel);

	//set pool sizes
	size_t locExpectedNumCombos = 100; //hopefully not often more than this
	dKinFitUtils->Set_MaxPoolSizes(Get_NumKinFitReactions(locEventLoop), locExpectedNumCombos);

	//pre-allocate matrix memory
	dKinFitUtils->Preallocate_MatrixMemory();

	return NOERROR;
}

size_t DKinFitResults_factory::Get_NumKinFitReactions(JEventLoop* locEventLoop)
{
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
		for(size_t loc_j = 0; loc_j < locReactionsSubset.size(); ++loc_j)
		{
			if(locReactionsSubset[loc_j]->Get_KinFitType() != d_NoFit)
				++locNumReactions;
		}
	}

	return locNumReactions;
}

//------------------
// evnt
//------------------
jerror_t DKinFitResults_factory::evnt(JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DKinFitResults_factory::evnt()");
#endif
	dConstraintResultsMap.clear();

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

	dKinFitter->Reset_NewEvent();
	for(size_t loc_i = 0; loc_i < locParticleCombos.size(); ++loc_i)
	{
		const DParticleCombo* locParticleCombo = locParticleCombos[loc_i];
		DKinFitType locKinFitType = locParticleCombo->Get_Reaction()->Get_KinFitType();
		if(locKinFitType == d_NoFit)
			continue; //don't do any kinematic fits!
		bool locP4IsFitFlag = ((locKinFitType != d_VertexFit) && (locKinFitType != d_SpacetimeFit));

		//Make DKinFitChain
		const DKinFitChain* locKinFitChain = dKinFitUtils->Make_KinFitChain(locParticleCombo, locKinFitType);

		//Make Constraints
		deque<DKinFitConstraint_Vertex*> locSortedVertexConstraints;
		set<DKinFitConstraint*> locConstraints = dKinFitUtils->Create_Constraints(locParticleCombo, locKinFitChain, locKinFitType, locSortedVertexConstraints);
		if(locConstraints.empty())
			continue; //Nothing to fit!

		//see if constraints (particles) are identical to a previous kinfit
		map<set<DKinFitConstraint*>, DKinFitResults*>::iterator locResultIterator = dConstraintResultsMap.find(locConstraints);
		if(locResultIterator != dConstraintResultsMap.end())
		{
			//this has been kinfit before, use the same result
			DKinFitResults* locKinFitResults = locResultIterator->second;
			if(locKinFitResults != NULL)
			{
				//previous kinfit succeeded, build the output DKinFitChain and register this combo with that fit
				set<DKinFitParticle*> locOutputKinFitParticles = locKinFitResults->Get_OutputKinFitParticles();
				const DKinFitChain* locOutputKinFitChain = dKinFitUtils->Build_OutputKinFitChain(locKinFitChain, locOutputKinFitParticles);
				locKinFitResults->Add_ParticleCombo(locParticleCombo, locOutputKinFitChain);
			}
			//else: the previous kinfit failed, so this one will too (don't save)
			continue;
		}

		//Set Constraint Initial Guesses //vertices, times, and missing p3's
		if(!locSortedVertexConstraints.empty())
			dKinFitUtils->Set_SpacetimeGuesses(locSortedVertexConstraints, locP4IsFitFlag);

		//Add constraints & perform fit
		dKinFitter->Reset_NewFit();
		dKinFitter->Add_Constraints(locConstraints);
		bool locFitStatus = dKinFitter->Fit_Reaction();

		//Build results (unless failed), and register
		DKinFitResults* locKinFitResults = NULL;
		if(locFitStatus)
		{
			set<DKinFitParticle*> locOutputKinFitParticles = dKinFitter->Get_KinFitParticles();
			const DKinFitChain* locOutputKinFitChain = dKinFitUtils->Build_OutputKinFitChain(locKinFitChain, locOutputKinFitParticles);
			locKinFitResults = Build_KinFitResults(locParticleCombo, locOutputKinFitChain);
		}
		dConstraintResultsMap[locConstraints] = locKinFitResults;
	}

	return NOERROR;
}

DKinFitResults* DKinFitResults_factory::Build_KinFitResults(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain)
{
	DKinFitResults* locKinFitResults = new DKinFitResults();
	locKinFitResults->Add_ParticleCombo(locParticleCombo, locKinFitChain);

	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	locKinFitResults->Set_KinFitType(locReaction->Get_KinFitType());

	locKinFitResults->Set_ConfidenceLevel(dKinFitter->Get_ConfidenceLevel());
	locKinFitResults->Set_ChiSq(dKinFitter->Get_ChiSq());
	locKinFitResults->Set_NDF(dKinFitter->Get_NDF());

	//locKinFitResults->Set_VEta(dKinFitter->Get_VEta());
	locKinFitResults->Set_VXi(dKinFitter->Get_VXi());
	//locKinFitResults->Set_V(dKinFitter->Get_V());

	locKinFitResults->Set_NumConstraints(dKinFitter->Get_NumConstraintEquations());
	locKinFitResults->Set_NumUnknowns(dKinFitter->Get_NumUnknowns());

	//Output particles and constraints
	locKinFitResults->Add_OutputKinFitParticles(dKinFitter->Get_KinFitParticles());
	locKinFitResults->Add_KinFitConstraints(dKinFitter->Get_KinFitConstraints());

	//Pulls

	//Build this:
	map<const JObject*, map<DKinFitPullType, double> > locPulls_JObject;

	//From this:
	map<DKinFitParticle*, map<DKinFitPullType, double> > locPulls_KinFitParticle;
	dKinFitter->Get_Pulls(locPulls_KinFitParticle);

	//By looping over the pulls:
	map<DKinFitParticle*, map<DKinFitPullType, double> >::iterator locMapIterator = locPulls_KinFitParticle.begin();
	for(; locMapIterator != locPulls_KinFitParticle.end(); ++locMapIterator)
	{
		DKinFitParticle* locOutputKinFitParticle = locMapIterator->first;
		DKinFitParticle* locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(locOutputKinFitParticle);
		const JObject* locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);

		locPulls_JObject[locSourceJObject] = locMapIterator->second;
	}
	//Set Pulls
	locKinFitResults->Set_Pulls(locPulls_JObject);

	//Particle Mapping
	//If any particles were NOT part of the kinematic fit, they are still added to the source -> output map
	set<DKinFitParticle*> locAllKinFitParticles = locKinFitChain->Get_AllParticles();
	set<DKinFitParticle*>::iterator locParticleIterator = locAllKinFitParticles.begin();
	for(; locParticleIterator != locAllKinFitParticles.end(); ++locParticleIterator)
	{
		const JObject* locSourceJObject = dKinFitUtils->Get_SourceJObject(*locParticleIterator);
		if(locSourceJObject != NULL)
		{
			locKinFitResults->Add_ParticleMapping_SourceToOutput(locSourceJObject, *locParticleIterator);
			continue; //*locParticleIterator was an input object //not directly used in the fit
		}
		DKinFitParticle* locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(*locParticleIterator);
		if(locInputKinFitParticle != NULL)
		{
			locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);
			if(locSourceJObject != NULL) //else was a decaying/missing particle: no source
				locKinFitResults->Add_ParticleMapping_SourceToOutput(locSourceJObject, *locParticleIterator);
		}
	}

	_data.push_back(locKinFitResults);
	return locKinFitResults;
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

