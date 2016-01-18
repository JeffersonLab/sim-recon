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
	size_t locExpectedNumCombos = 50;
	dKinFitUtils->Set_MaxPoolSizes(Get_NumReactions(locEventLoop), locExpectedNumCombos);

	//pre-allocate matrix memory
	dKinFitUtils->Preallocate_MatrixMemory();

	return NOERROR;
}

size_t DKinFitResults_factory::Get_NumReactions(JEventLoop* locEventLoop)
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
		locNumReactions += locReactionsSubset.size();
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

		//Make DKinFitChain
		const DKinFitChain* locKinFitChain = dKinFitUtils->Make_KinFitChain(locParticleCombo, locKinFitType);

		//Make Constraints
		deque<DKinFitConstraint_Vertex*> locSortedVertexConstraints;
		set<DKinFitConstraint*> locConstraints = dKinFitUtils->Create_Constraints(locKinFitChain, locKinFitType, locSortedVertexConstraints);

		//see if constraints (particles) are identical to a previous kinfit
		map<set<DKinFitConstraint*>, DKinFitResults*>::iterator locResultIterator = dConstraintResultsMap.find(locConstraints);
		if(locResultIterator != dConstraintResultsMap.end())
		{
			//this has been kinfit before, use the same result
			DKinFitResults* locKinFitResults = locResultIterator->second;
			if(locKinFitResults != NULL) //if false: the previous kinfit failed, this one will too
				locKinFitResults->Add_ParticleCombo(locParticleCombo, locKinFitChain); //previous kinfit succeeded, register this combo with that fit
			continue;
		}

		//Set Constraint Initial Guesses //vertices, times, and missing p3's
		if(!locSortedVertexConstraints.empty())
			dKinFitUtils->Set_SpacetimeGuesses(locSortedVertexConstraints);

		//Add constraints & perform fit
		dKinFitter->Reset_NewFit();
		dKinFitter->Add_Constraints(locConstraints);
		bool locFitStatus = dKinFitter->Fit_Reaction();

		//Build results (unless failed), and register
		DKinFitResults* locKinFitResults = locFitStatus ? Build_KinFitResults(locParticleCombo, locKinFitChain) : NULL;
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

	locKinFitResults->Set_VEta(dKinFitter->Get_VEta());
	locKinFitResults->Set_VXi(dKinFitter->Get_VXi());
	locKinFitResults->Set_V(dKinFitter->Get_V());

	locKinFitResults->Set_NumConstraints(dKinFitter->Get_NumConstraintEquations());
	locKinFitResults->Set_NumUnknowns(dKinFitter->Get_NumUnknowns());

	//Output particles and constraints
	locKinFitResults->Add_OutputKinFitParticles(dKinFitter->Get_KinFitParticles());
	locKinFitResults->Add_KinFitConstraints(dKinFitter->Get_KinFitConstraints());

	//Particle Maps & Pulls
	//Build this:
	map<const JObject*, map<DKinFitPullType, double> > locPulls_JObject;
	//From this:
	map<DKinFitParticle*, map<DKinFitPullType, double> > locPulls_KinFitParticle;
	dKinFitter->Get_Pulls(locPulls_KinFitParticle);
	//By looping over the pulls (and can do mapping too): //only particles with pulls need mapping (others have no source objects)
	map<DKinFitParticle*, map<DKinFitPullType, double> >::iterator locMapIterator = locPulls_KinFitParticle.begin();
	for(; locMapIterator != locPulls_KinFitParticle.end(); ++locMapIterator)
	{
		DKinFitParticle* locOutputKinFitParticle = locMapIterator->first;
		DKinFitParticle* locInputKinFitParticle = dKinFitUtils->Get_InputKinFitParticle(locOutputKinFitParticle);
		const JObject* locSourceJObject = dKinFitUtils->Get_SourceJObject(locInputKinFitParticle);

		locPulls_JObject[locSourceJObject] = locMapIterator->second;
		locKinFitResults->Add_ParticleMapping_SourceToInput(locSourceJObject, locInputKinFitParticle);
		locKinFitResults->Add_ParticleMapping_InputToOutput(locInputKinFitParticle, locOutputKinFitParticle);
	}
	//Set Pulls
	locKinFitResults->Set_Pulls(locPulls_JObject);

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

