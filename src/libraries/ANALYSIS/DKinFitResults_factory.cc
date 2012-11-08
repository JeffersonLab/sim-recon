#include "DKinFitResults_factory.h"

//------------------
// init
//------------------
jerror_t DKinFitResults_factory::init(void)
{
	dDebugLevel = 0;
	dKinFitDebugLevel = 0;
	dLinkVerticesFlag = false;
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
	dKinFitter.Set_BField(locApplication->GetBfield());

	gPARMS->SetDefaultParameter("KINFIT:KINFIT_DEBUGLEVEL", dKinFitDebugLevel);
	gPARMS->SetDefaultParameter("KINFIT:DEBUGLEVEL", dDebugLevel);
	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);

	dKinFitter.Set_DebugLevel(dKinFitDebugLevel);
	dKinFitter.Set_LinkVerticesFlag(dLinkVerticesFlag);

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

		//check previous particle combos for duplicates: if all of the steps are the same AND the fit type is the same, then the kinfitresults will be the same: just copy!
			//this would be used if you want to perform the exact same kinematic fit, but want to perform different DAnalyssisActions on them
		bool locComboMatchFlag = false;
		for(size_t loc_j = 0; loc_j < _data.size(); ++loc_j)
		{
			const DParticleCombo* locPreviousParticleCombo = _data[loc_j]->Get_ParticleCombo();
			if(!locParticleCombo->Will_KinFitBeIdentical(locPreviousParticleCombo))
				continue;
			if(dDebugLevel > 0)
				cout << "kinfit results will be the same: clone and save" << endl;
			//kinfit results will be the same: clone and save
			locComboMatchFlag = true;
			DKinFitResults* locKinFitResults = new DKinFitResults(*(_data[loc_j]));
			locKinFitResults->Set_ParticleCombo(locParticleCombo);
			_data.push_back(locKinFitResults);
			break;
		}
		if(locComboMatchFlag)
			continue; //kinfit results already saved, continue;

		deque<deque<const DKinFitParticle*> > locInitialKinFitParticles, locFinalKinFitParticles;
		Setup_KinFit(locParticleCombo, locInitialKinFitParticles, locFinalKinFitParticles);
		if(dDebugLevel > 0)
			cout << "Perform Primary Kinematic Fit" << endl;
		if(dKinFitter.Fit_Reaction()) //if fit fails: no kinfit results, BUT will still generate new DParticleCombo (using old info though!)
			Build_KinFitResults(locParticleCombo, locInitialKinFitParticles, locFinalKinFitParticles);
	}

	return NOERROR;
}

void DKinFitResults_factory::Setup_KinFit(const DParticleCombo* locParticleCombo, deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles)
{
	if(dDebugLevel > 0)
		cout << "Setup Primary Kinematic Fit" << endl;

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
				locFinalKinFitParticles_Step.push_back(locKinFitParticle);
				locDecayingParticleStepMap[locDecayStepIndex] = locKinFitParticle;
			}
			else if(locKinematicData->charge() == 0) //detected neutral shower
			{
				const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);
				locKinFitParticle = dKinFitter.Make_DetectedParticle(locNeutralParticleHypothesis, (locKinFitType == d_P4Fit));
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

	bool locUseRFTimeFlag;
	double locRFTime = 0.0;
	double locRFUncertainty = 0.0;

	//VERTEX OR SPACETIME: Group particles by detached vertex (one deque for each constraint/vertex), calculate initial guesses for vertices & times
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_Vertices, locFinalKinFitParticles_Vertices;
	deque<TVector3> locVertexGuesses; //guesses for each vertex: kinfit (input is crude guess) (if kinfit fails, use crude guess)
	deque<double> locTimeGuesses;
	deque<bool> locUseRFTimeFlags;
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Setup Vertex & Spacetime Constraints" << endl;
		deque<size_t> locStepIndicesToHandle_Vertex;
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
			locStepIndicesToHandle_Vertex.push_back(loc_i);
		while(!locStepIndicesToHandle_Vertex.empty())
		{
			deque<const DKinFitParticle*> locInitialKinFitParticles_Vertex;
			deque<const DKinFitParticle*> locFinalKinFitParticles_Vertex;
			deque<size_t> locIncludedStepIndices;
			size_t locStepIndex = locStepIndicesToHandle_Vertex[0];

			Setup_VertexConstraint(locParticleCombo, locStepIndicesToHandle_Vertex[0], locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, locIncludedStepIndices);

			//remove steps included in the vertex constraint from the to-handle deque
			for(size_t loc_i = 0; loc_i < locIncludedStepIndices.size(); ++loc_i)
			{
				for(deque<size_t>::iterator locIterator = locStepIndicesToHandle_Vertex.begin(); locIterator != locStepIndicesToHandle_Vertex.end(); ++locIterator)
				{
					if(locIncludedStepIndices[loc_i] == (*locIterator))
					{
						locStepIndicesToHandle_Vertex.erase(locIterator);
						break;
					}
				}
			}

			//make sure there's enough detected particles to constrain the vertex //undeterminable here if dLinkVerticesFlag = true
				//may be ok with 0 or 1 tracks if enough decaying particles and they're constrained elsewhere (will figure it out later)
			size_t locNumDetectedChargedParticles = 0;
			for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_Vertex.size(); ++loc_i)
			{
				if((locFinalKinFitParticles_Vertex[loc_i]->Get_KinFitParticleType() == d_DetectedParticle) && (locFinalKinFitParticles_Vertex[loc_i]->Get_Charge() != 0))
					++locNumDetectedChargedParticles;
			}
			if((locNumDetectedChargedParticles < 2) && (!dLinkVerticesFlag))
				continue; //not enough detected tracks!!!

			locInitialKinFitParticles_Vertices.push_back(locInitialKinFitParticles_Vertex);
			locFinalKinFitParticles_Vertices.push_back(locFinalKinFitParticles_Vertex);

			TVector3 locVertexGuess = Calc_VertexGuess(locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex);
			locVertexGuesses.push_back(locVertexGuess);
			if((locKinFitType == d_VertexFit) || (locKinFitType == d_P4AndVertexFit)) //vertex
				continue;

			//spacetime
			locUseRFTimeFlag = ((locStepIndex == 0) && locFirstParticleIsBeamFlag);
			locUseRFTimeFlags.push_back(locUseRFTimeFlag);
			double locTimeGuess = Calc_TimeGuess(locFinalKinFitParticles_Vertex, DVector3(locVertexGuess.X(),locVertexGuess.Y(),locVertexGuess.Z()), locUseRFTimeFlag, locRFTime);
			locTimeGuesses.push_back(locTimeGuess);
		}
	}


	//ADD CONSTRAINTS
	dKinFitter.Reset_NewFit();

	//P4: Create and add constraints
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Setup P4 Constraints" << endl;
		deque<size_t> locStepIndicesToHandle_P4;
		for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
			locStepIndicesToHandle_P4.push_back(loc_i);
		while(!locStepIndicesToHandle_P4.empty())
		{
			deque<const DKinFitParticle*> locInitialKinFitParticles_P4;
			deque<const DKinFitParticle*> locFinalKinFitParticles_P4;
			deque<size_t> locIncludedStepIndices;

			Setup_P4Constraint(locParticleCombo, locStepIndicesToHandle_P4[0], locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locIncludedStepIndices);

			//remove steps included in the p4 constraint from the to-handle deque
			for(size_t loc_i = 0; loc_i < locIncludedStepIndices.size(); ++loc_i)
			{
				for(deque<size_t>::iterator locIterator = locStepIndicesToHandle_P4.begin(); locIterator != locStepIndicesToHandle_P4.end(); ++locIterator)
				{
					if(locIncludedStepIndices[loc_i] == (*locIterator))
					{
						locStepIndicesToHandle_P4.erase(locIterator);
						break;
					}
				}
			}

			if(dDebugLevel > 0)
			{
				cout << "INIT particles in p4 constraint: " << endl;
				for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_P4.size(); ++loc_i)
					cout << locInitialKinFitParticles_P4[loc_i]->Get_Mass() << ", " << locInitialKinFitParticles_P4[loc_i]->Get_Charge() << endl;
				cout << "FINAL particles in p4 constraint: " << endl;
				for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_P4.size(); ++loc_i)
					cout << locFinalKinFitParticles_P4[loc_i]->Get_Mass() << ", " << locFinalKinFitParticles_P4[loc_i]->Get_Charge() << endl;
			}
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Create P4 Constraint" << endl;
			dKinFitter.Add_P4Constraint(locInitialKinFitParticles_P4, locFinalKinFitParticles_P4);
		}
	}

	//VERTEX: Add constraints
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_P4AndVertexFit))
	{
		Remove_BadVertexConstraints(locInitialKinFitParticles_Vertices, locFinalKinFitParticles_Vertices);
		for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_Vertices.size(); ++loc_i)
			dKinFitter.Add_VertexConstraint(locInitialKinFitParticles_Vertices[loc_i], locFinalKinFitParticles_Vertices[loc_i], locVertexGuesses[loc_i]);
	}

	//SPACETIME: Add constraints, and RF time if needed
	if((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_Vertices.size(); ++loc_i)
			dKinFitter.Add_SpacetimeConstraint(locInitialKinFitParticles_Vertices[loc_i], locFinalKinFitParticles_Vertices[loc_i], locUseRFTimeFlags[loc_i], locVertexGuesses[loc_i], locTimeGuesses[loc_i]);

		if(locFirstParticleIsBeamFlag)
			dKinFitter.Set_RFTime(locRFTime, locRFUncertainty, locInitialKinFitParticles[0][0]);
	}
}

void DKinFitResults_factory::Remove_BadVertexConstraints(deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Vertices, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Vertices) const
{
	//Exclude vertex fits that don't have enough particles
		//this is assumed to already have been checked if dLinkVerticesFlag = false (trivial), but not if true
			//this is because could have a decaying particle defined in one fit and used to constrain another

	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_OKVertices, locFinalKinFitParticles_OKVertices;
	if(!dLinkVerticesFlag)
		return; //already done prior to here!!!
//THIS ROUTINE IS UNTESTED!!!
	const DKinFitParticle* locKinFitParticle;
	deque<const DKinFitParticle*> locConstrainedDecayingParticles;
	while(!locFinalKinFitParticles_Vertices.empty())
	{
		//loop over the vertex fits, saving the good ones as "OK"
		deque<size_t> locSavedVertexFits;
		for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_Vertices.size(); ++loc_i)
		{
			//find out how many of the final particles in this fit are constrained
			size_t locNumConstrainedParticles = 0;
			for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_Vertices[loc_i].size(); ++loc_j)
			{
				locKinFitParticle = locFinalKinFitParticles_Vertices[loc_i][loc_j];
				if((locKinFitParticle->Get_KinFitParticleType() == d_DetectedParticle) && (locKinFitParticle->Get_Charge() != 0))
					++locNumConstrainedParticles; //detected charged particle
				else if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle)
				{
					for(size_t loc_k = 0; loc_k < locConstrainedDecayingParticles.size(); ++loc_k)
					{
						if(locConstrainedDecayingParticles[loc_k] == locKinFitParticle)
						{
							++locNumConstrainedParticles; //decaying particle constrained by another constraint
							break;
						}
					}
				}
			}
			//find out if the initial particle in this fit is constrained
			if(!locInitialKinFitParticles_Vertices[loc_i].empty())
			{
				locKinFitParticle = locInitialKinFitParticles_Vertices[loc_i][0];
				if(locKinFitParticle != NULL)
				{
					if(locKinFitParticle->Get_KinFitParticleType() == d_BeamParticle) //beam photon
					{
						if((*locKinFitParticle->Get_CovarianceMatrix())(3, 3) > 0.0) //constrained!!
							++locNumConstrainedParticles; //beam particle has non-zero uncertainty on beam position
					}
					else if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle) //decaying particle
					{
						for(size_t loc_k = 0; loc_k < locConstrainedDecayingParticles.size(); ++loc_k)
						{
							if(locConstrainedDecayingParticles[loc_k] == locKinFitParticle)
							{
								++locNumConstrainedParticles; //decaying particle constrained by another constraint
								break;
							}
						}
					}
				}
			}

			//if enough constrained particles, this fit is OK
			if(locNumConstrainedParticles >= 2)
			{
				//label any other decaying particles in this fit as constrained
				for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_Vertices[loc_i].size(); ++loc_j)
				{
					locKinFitParticle = locFinalKinFitParticles_Vertices[loc_i][loc_j];
					if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
						continue;
					bool locParticleConstrainedFlag = false;
					for(size_t loc_k = 0; loc_k < locConstrainedDecayingParticles.size(); ++loc_k)
					{
						if(locConstrainedDecayingParticles[loc_k] == locKinFitParticle)
						{
							locParticleConstrainedFlag = true;
							break;
						}
					}
					if(!locParticleConstrainedFlag)
						locConstrainedDecayingParticles.push_back(locKinFitParticle);
				}
				//label the decaying initial particle as constrained (if it exists)
				if(!locInitialKinFitParticles_OKVertices[loc_i].empty())
				{
					locKinFitParticle = locInitialKinFitParticles_Vertices[loc_i][0];
					if(locKinFitParticle != NULL)
					{
						if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle) //decaying particle
						{
							bool locParticleConstrainedFlag = false;
							for(size_t loc_k = 0; loc_k < locConstrainedDecayingParticles.size(); ++loc_k)
							{
								if(locConstrainedDecayingParticles[loc_k] == locKinFitParticle)
								{
									locParticleConstrainedFlag = true;
									break;
								}
							}
							if(!locParticleConstrainedFlag)
								locConstrainedDecayingParticles.push_back(locKinFitParticle);
						}
					}
				}

				//mark vertex fit as OK
				locInitialKinFitParticles_OKVertices.push_back(locInitialKinFitParticles_Vertices[loc_i]);
				locFinalKinFitParticles_OKVertices.push_back(locFinalKinFitParticles_Vertices[loc_i]);
				locSavedVertexFits.push_back(loc_i);
			}
		}
		if(locSavedVertexFits.empty())
			break; //no progress: all remaining fits are bad

		for(int loc_i = locSavedVertexFits.size() - 1; loc_i >= 0; --loc_i)
		{
			locInitialKinFitParticles_Vertices.erase(locInitialKinFitParticles_Vertices.begin() + locSavedVertexFits[loc_i]);
			locFinalKinFitParticles_Vertices.erase(locFinalKinFitParticles_Vertices.begin() + locSavedVertexFits[loc_i]);
		}
	}
	locInitialKinFitParticles_Vertices = locInitialKinFitParticles_OKVertices;
	locFinalKinFitParticles_Vertices = locFinalKinFitParticles_OKVertices;
}

void DKinFitResults_factory::Setup_P4Constraint(const DParticleCombo* locParticleCombo, size_t locStepIndex, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles, deque<const DKinFitParticle*>& locInitialKinFitParticles_P4, deque<const DKinFitParticle*>& locFinalKinFitParticles_P4, deque<size_t>& locIncludedStepIndices)
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
		if(!locParticleCombo->Check_IfDecayingParticleExcludedFromP4KinFit(locStepIndex))
		{
			if(IsFixedMass(locPID)) //else is decaying resonance particle: don't constrain (don't recurse either (already done!))
				locInitialKinFitParticles_P4.push_back(locInitialKinFitParticles[locStepIndex][0]);
		}
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
		{
			if(IsFixedMass(locPID))
				locFinalKinFitParticles_P4.push_back(locFinalKinFitParticles[locStepIndex][loc_j]);
		}
		else if(locDecayStepIndex >= 0) //decaying particle
		{
			if(locParticleCombo->Check_IfDecayingParticleExcludedFromP4KinFit(locDecayStepIndex)) //go to the next step!!
				Setup_P4Constraint(locParticleCombo, locDecayStepIndex, locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locIncludedStepIndices);
			else if(!IsFixedMass(locPID)) //go to the next step!!
				Setup_P4Constraint(locParticleCombo, locDecayStepIndex, locInitialKinFitParticles, locFinalKinFitParticles, locInitialKinFitParticles_P4, locFinalKinFitParticles_P4, locIncludedStepIndices);
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

TVector3 DKinFitResults_factory::Calc_VertexGuess(const deque<const DKinFitParticle*>& locInitialKinFitParticles, const deque<const DKinFitParticle*>& locFinalKinFitParticles)
{
	//choose midpoint of DOCA line as init guess for vertex
	//do vertex-only kinfit to find the super-init kinfit guess
	deque<const DKinFitParticle*> locVertexFindParticles;
	DKinFitParticleType locKinFitParticleType;
	for(size_t loc_i = 0; loc_i < locFinalKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = locFinalKinFitParticles[loc_i]->Get_KinFitParticleType();
		if((locKinFitParticleType == d_DetectedParticle) && (locFinalKinFitParticles[loc_i]->Get_Charge() != 0))
			locVertexFindParticles.push_back(locFinalKinFitParticles[loc_i]);
	}
	DVector3 locTempInitVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexFindParticles);
	TVector3 locInitVertex(locTempInitVertex.X(),locTempInitVertex.Y(),locTempInitVertex.Z());

	//make sure there's enough detected particles to constrain the vertex for the guess
	size_t locNumDetectedChargedParticles = 0;
	for(size_t loc_i = 0; loc_i < locFinalKinFitParticles.size(); ++loc_i)
	{
		if((locFinalKinFitParticles[loc_i]->Get_KinFitParticleType() == d_DetectedParticle) && (locFinalKinFitParticles[loc_i]->Get_Charge() != 0))
			++locNumDetectedChargedParticles;
	}
	if(locNumDetectedChargedParticles < 2)
	  return locInitVertex; //not enough detected tracks to fit

	dKinFitter.Reset_NewFit();
	const DKinFitConstraint_Vertex* locVertexConstraint = dKinFitter.Add_VertexConstraint(locInitialKinFitParticles, locFinalKinFitParticles, locInitVertex);

	if(dDebugLevel > 0)
		cout << "Perform init vertex guess kinematic fit" << endl;
	bool locFitStatus = dKinFitter.Fit_Reaction();
	if(dDebugLevel > 0)
		cout << "Init vertex guess kinematic fit finished." << endl;
	if(locFitStatus)
		return locVertexConstraint->Get_CommonVertex();
	return locInitVertex;
//COULD PROPAGATE TRACK INFO TO POCA TO VERTEX RIGHT HERE!!! //advantage over using the kinfitter to do it: take energy loss into account
}

double DKinFitResults_factory::Calc_TimeGuess(const deque<const DKinFitParticle*>& locFinalKinFitParticles, DVector3 locVertexGuess, bool locUseRFTimeFlag, double locRFTime)
{
	//if RF: propagate RF time to the vertex-z of the init vertex guess
	//if no RF:
		//propagate each track time to the DOCA to the init vertex guess
		//take avg of all track times, weighted by their uncertainties

	if(locUseRFTimeFlag)
	{
		DVector3 locRFPosition(0.0, 0.0, 65.0); //i know this is bad... but really should have a DRFBunch object...
		return (locRFTime + (locVertexGuess.Z() - locRFPosition.Z())/29.9792458);
	}

	deque<const DKinFitParticle*> locTimeFindParticles;
	DKinFitParticleType locKinFitParticleType;
	for(size_t loc_i = 0; loc_i < locFinalKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = locFinalKinFitParticles[loc_i]->Get_KinFitParticleType();
		if((locKinFitParticleType == d_DetectedParticle) && (locFinalKinFitParticles[loc_i]->Get_Charge() != 0))
			locTimeFindParticles.push_back(locFinalKinFitParticles[loc_i]);
	}

	return dAnalysisUtilities->Calc_CrudeTime(locTimeFindParticles, locVertexGuess);
}

void DKinFitResults_factory::Build_KinFitResults(const DParticleCombo* locParticleCombo, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Input, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Input)
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	DKinFitResults* locKinFitResults = new DKinFitResults();
	locKinFitResults->Set_ParticleCombo(locParticleCombo);

	locKinFitResults->Set_KinFitName(locReaction->Get_ReactionName());
	locKinFitResults->Set_KinFitType(locReaction->Get_KinFitType());

	locKinFitResults->Set_ConfidenceLevel(dKinFitter.Get_ConfidenceLevel());
	locKinFitResults->Set_ChiSq(dKinFitter.Get_ChiSq());
	locKinFitResults->Set_NDF(dKinFitter.Get_NDF());
	locKinFitResults->Set_VEta(dKinFitter.Get_VEta());
	locKinFitResults->Set_VXi(dKinFitter.Get_VXi());

	//save the constraints
	deque<const DKinFitConstraint*> locKinFitConstraints;
	dKinFitter.Get_KinFitConstraints(locKinFitConstraints);
	locKinFitResults->Set_KinFitConstraints(locKinFitConstraints);

	//create output deques from input deques & mapping
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_Output, locFinalKinFitParticles_Output;
	locInitialKinFitParticles_Output.resize(locInitialKinFitParticles_Input.size());
	const DKinFitParticle* locOutputKinFitParticle;
	for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_Input.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locInitialKinFitParticles_Input[loc_i].size(); ++loc_j)
		{
			locOutputKinFitParticle = dKinFitter.Get_OutputKinFitParticle(locInitialKinFitParticles_Input[loc_i][loc_j]);
			locInitialKinFitParticles_Output[loc_i].push_back(locOutputKinFitParticle);
		}
	}
	locFinalKinFitParticles_Output.resize(locFinalKinFitParticles_Input.size());
	for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_Input.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_Input[loc_i].size(); ++loc_j)
		{
			locOutputKinFitParticle = dKinFitter.Get_OutputKinFitParticle(locFinalKinFitParticles_Input[loc_i][loc_j]);
			locFinalKinFitParticles_Output[loc_i].push_back(locOutputKinFitParticle);
		}
	}
	locKinFitResults->Set_InitialKinFitParticles(locInitialKinFitParticles_Output);
	locKinFitResults->Set_FinalKinFitParticles(locFinalKinFitParticles_Output);

	//convert the map of the data: convert from (input particle -> data) to (output particle -> data)
	map<const DKinFitParticle*, const DKinematicData*> locParticleMapping_Output;
	dKinFitter.Get_ParticleMapping_OutputToSource(locParticleMapping_Output);
	locKinFitResults->Set_ParticleMapping(locParticleMapping_Output);

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

