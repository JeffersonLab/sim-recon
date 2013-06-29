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
	dKinFitter.Set_BField(locApplication->GetBfield());

	gPARMS->SetDefaultParameter("KINFIT:KINFITDEBUGLEVEL", dKinFitDebugLevel);
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
			//this would be used if you want to perform the exact same kinematic fit, but want to perform different DAnalysisActions on them
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

		map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > > locDecayingKinFitParticles;
		if(!Setup_KinFit(locEventLoop, locParticleCombo, locDecayingKinFitParticles))
		{
			if(dDebugLevel > 0)
				cout << "Kinematic Fit Setup Failed." << endl;
			continue;
		}
		if(dDebugLevel > 0)
			cout << "Perform Primary Kinematic Fit" << endl;
		if(dKinFitter.Fit_Reaction()) //if fit fails: no kinfit results, BUT will still generate new DParticleCombo (using old info though!)
			Build_KinFitResults(locParticleCombo, locDecayingKinFitParticles);
	}

	return NOERROR;
}

bool DKinFitResults_factory::Setup_KinFit(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo, map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locDecayingKinFitParticles)
{
	if(dDebugLevel > 0)
		cout << "Setup Primary Kinematic Fit" << endl;

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
	deque<const DKinFitParticle*> locUndetectedParticles_NoP4Guess;
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
			if(IsFixedMass(locPID))
				locUndetectedParticles_NoP4Guess.push_back(locKinFitParticle);
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
				if(IsFixedMass(locPID))
					locUndetectedParticles_NoP4Guess.push_back(locKinFitParticle);
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
				if(IsFixedMass(locPID))
					locUndetectedParticles_NoP4Guess.push_back(locKinFitParticle);
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

	//P4: Setup constraint info
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_P4s, locFinalKinFitParticles_P4s;
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
				cout << "INIT particles in p4 constraint (mass, q): " << endl;
				for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_P4.size(); ++loc_i)
					cout << locInitialKinFitParticles_P4[loc_i]->Get_Mass() << ", " << locInitialKinFitParticles_P4[loc_i]->Get_Charge() << endl;
				cout << "FINAL particles in p4 constraint (mass, q): " << endl;
				for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_P4.size(); ++loc_i)
					cout << locFinalKinFitParticles_P4[loc_i]->Get_Mass() << ", " << locFinalKinFitParticles_P4[loc_i]->Get_Charge() << endl;
			}
			locInitialKinFitParticles_P4s.push_back(locInitialKinFitParticles_P4);
			locFinalKinFitParticles_P4s.push_back(locFinalKinFitParticles_P4);
		}
	}

	//GET P4 GUESSES OF DECAYING PARTICLES (for use in guessing initial vertex positions)
	//loop through the p4 constraint info, find constraints which contain only one missing or decaying particle: used as starting point for assigning them to different constraints
	map<const DKinFitParticle*, TVector3> locP4GuessMap;
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		deque<pair<const DKinFitParticle*, size_t> > locConstrainableParticles; //size_t is constraint index (index of locFinalKinFitParticles_P4s)
		deque<const DKinFitParticle*> locConstrainedParticles;
		deque<size_t> locConstraintsSetIndices;
		while(locConstrainedParticles.size() < locUndetectedParticles_NoP4Guess.size())
		{
			if(locConstrainableParticles.empty())
			{
				if(!Find_ConstrainableParticles(locConstrainableParticles, locConstrainedParticles, locInitialKinFitParticles_P4s, locFinalKinFitParticles_P4s, locConstraintsSetIndices))
					return false; //more decaying/missing particles than constraints, or circular dependency: cannot fit
			}

			while(!locConstrainableParticles.empty())
			{
				bool locAlreadyConstrainedFlag = false;
				//see if the particle is already constrained (e.g. constrainable in more than one constraint)
				for(size_t loc_j = 0; loc_j < locConstrainedParticles.size(); ++loc_j)
				{
					if(locConstrainableParticles.back().first != locConstrainedParticles[loc_j])
						continue;
					locAlreadyConstrainedFlag = true;
					break;
				}
				if(locAlreadyConstrainedFlag)
				{
					locConstrainableParticles.pop_back();
					continue;
				}

				TVector3 locMomentum;
				Calc_P4Guess(locConstrainableParticles.back(), locInitialKinFitParticles_P4s, locFinalKinFitParticles_P4s, locMomentum);
				locP4GuessMap[locConstrainableParticles.back().first] = locMomentum;
				locConstraintsSetIndices.push_back(locConstrainableParticles.back().second);
				locConstrainedParticles.push_back(locConstrainableParticles.back().first);
				locConstrainableParticles.pop_back();
			}
		}
	}


	//VERTEX OR SPACETIME: Group particles by detached vertex (one deque for each constraint/vertex)
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_Vertices_Staging, locFinalKinFitParticles_Vertices_Staging;
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

			//make sure there's enough detected particles to constrain the vertex //not-completely-determinable here if dLinkVerticesFlag = true
				//may be ok with 0 or 1 tracks if enough decaying particles and they're constrained elsewhere (will figure it out later)
			size_t locNumDetectedChargedParticles = 0;
			for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_Vertex.size(); ++loc_i)
			{
				if((locFinalKinFitParticles_Vertex[loc_i]->Get_KinFitParticleType() == d_DetectedParticle) && (locFinalKinFitParticles_Vertex[loc_i]->Get_Charge() != 0))
					++locNumDetectedChargedParticles;
			}
			if((locNumDetectedChargedParticles < 2) && (!dLinkVerticesFlag))
				continue; //not enough detected tracks!!!

			locInitialKinFitParticles_Vertices_Staging.push_back(locInitialKinFitParticles_Vertex);
			locFinalKinFitParticles_Vertices_Staging.push_back(locFinalKinFitParticles_Vertex);
		}
	}

	//Calculate initial guesses for vertices & times
		//if >= 2 charged tracks then do kinfit
		//if < 2 charged tracks but decaying track, do a kinfit with 
	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_Vertices, locFinalKinFitParticles_Vertices;
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitResults_factory: Calc Vertex & Time Guesses" << endl;
		int locVertexFindFlag = 0; //1 after first guess
		map<const DKinFitParticle*, TVector3> locDecayingParticleTrackPointGuesses;
		size_t locNumStagedConstraints = 0;
		do //loop of loops: first calc all vertex guesses for >= 2 charged tracks, then use this position information for decaying particles in fits with < 2 charged tracks
		{
			locNumStagedConstraints = locFinalKinFitParticles_Vertices_Staging.size();
			deque<deque<const DKinFitParticle*> >::iterator locIterator_Initial = locInitialKinFitParticles_Vertices_Staging.begin();
			deque<deque<const DKinFitParticle*> >::iterator locIterator_Final = locFinalKinFitParticles_Vertices_Staging.begin();
			while(locIterator_Final != locFinalKinFitParticles_Vertices_Staging.end()) //loop over all vertex constraints still staged
			{
				deque<const DKinFitParticle*> locInitialKinFitParticles_Vertex = *locIterator_Initial;
				deque<const DKinFitParticle*> locFinalKinFitParticles_Vertex = *locIterator_Final;
				TVector3 locVertexGuess;
				deque<const DKinFitParticle*> locFinalKinFitParticles_Vertex_Updated; //propagate the track information if charged & detected
				bool locVertexOKFlag = Calc_VertexGuess(locEventLoop, locInitialKinFitParticles_Vertex, locFinalKinFitParticles_Vertex, locVertexGuess, locFinalKinFitParticles_Vertex_Updated, locVertexFindFlag, locP4GuessMap, locDecayingParticleTrackPointGuesses);

				if(!locVertexOKFlag)
				{
					++locIterator_Final;
					++locIterator_Initial;
					continue;
				}

				//update the master kinfitparticle deques to use the new kinfitparticles (track info propagated)
				for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_Vertex_Updated.size(); ++loc_i)
				{
					if(locFinalKinFitParticles_Vertex_Updated[loc_i] == NULL)
						continue; //not updated
					//else updated:
					for(size_t loc_j = 0; loc_j < locFinalKinFitParticles.size(); ++loc_j)
					{
						bool locMatchFoundFlag = false;
						for(size_t loc_k = 0; loc_k < locFinalKinFitParticles[loc_j].size(); ++loc_k)
						{
							if(locFinalKinFitParticles_Vertex[loc_i] != locFinalKinFitParticles[loc_j][loc_k])
								continue;
							locFinalKinFitParticles[loc_j][loc_k] = locFinalKinFitParticles_Vertex_Updated[loc_i]; //update master deque
							locFinalKinFitParticles_Vertex[loc_i] = locFinalKinFitParticles_Vertex_Updated[loc_i];
							locMatchFoundFlag = true;
							break;
						}
						if(locMatchFoundFlag)
							break;
					}
				}

				//guess obtained: unstage this constraint
				locIterator_Final = locFinalKinFitParticles_Vertices_Staging.erase(locIterator_Final);
				locIterator_Initial = locInitialKinFitParticles_Vertices_Staging.erase(locIterator_Initial);

				//save constraint info
				locInitialKinFitParticles_Vertices.push_back(locInitialKinFitParticles_Vertex);
				locFinalKinFitParticles_Vertices.push_back(locFinalKinFitParticles_Vertex);

				locVertexGuesses.push_back(locVertexGuess);
				if((locKinFitType == d_VertexFit) || (locKinFitType == d_P4AndVertexFit)) //vertex
					continue;

				//spacetime
				locUseRFTimeFlag = false;
				if(locFirstParticleIsBeamFlag && (!locInitialKinFitParticles_Vertex.empty()))
					locUseRFTimeFlag = ((locInitialKinFitParticles_Vertex[0]->Get_Charge() == 0) && (!(locInitialKinFitParticles_Vertex[0]->Get_Mass() > 0.0))); //true if photon
				locUseRFTimeFlags.push_back(locUseRFTimeFlag);
				double locTimeGuess = Calc_TimeGuess(locFinalKinFitParticles_Vertex, DVector3(locVertexGuess.X(), locVertexGuess.Y(), locVertexGuess.Z()), locUseRFTimeFlag, locRFTime);
				locTimeGuesses.push_back(locTimeGuess);
			}
			locVertexFindFlag = 1; //now use decaying particle info (if present)
		}
		while(locFinalKinFitParticles_Vertices_Staging.size() < locNumStagedConstraints); //if no progress then break
	}

	//ADD CONSTRAINTS
	dKinFitter.Reset_NewFit();

	//P4: ADD CONSTRAINTS
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		for(size_t loc_i = 0; loc_i < locFinalKinFitParticles_P4s.size(); ++loc_i)
		{
			if(dDebugLevel > 10)
				cout << "DKinFitResults_factory: Create P4 Constraint" << endl;
			dKinFitter.Add_P4Constraint(locInitialKinFitParticles_P4s[loc_i], locFinalKinFitParticles_P4s[loc_i]);
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
		Remove_BadVertexConstraints(locInitialKinFitParticles_Vertices, locFinalKinFitParticles_Vertices);
		for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_Vertices.size(); ++loc_i)
			dKinFitter.Add_SpacetimeConstraint(locInitialKinFitParticles_Vertices[loc_i], locFinalKinFitParticles_Vertices[loc_i], locUseRFTimeFlags[loc_i], locVertexGuesses[loc_i], locTimeGuesses[loc_i]);
		if(locFirstParticleIsBeamFlag)
			dKinFitter.Set_RFTime(locRFTime, locRFUncertainty, locInitialKinFitParticles[0][0]);
	}

	return true;
}

void DKinFitResults_factory::Remove_BadVertexConstraints(deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_Vertices, deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_Vertices) const
{
	//Exclude vertex fits that don't have enough particles
		//this is assumed to already have been checked if dLinkVerticesFlag = false (trivial), but not if true
			//this is because could have a decaying particle defined in one fit and used to constrain another

	deque<deque<const DKinFitParticle*> > locInitialKinFitParticles_OKVertices, locFinalKinFitParticles_OKVertices;
	if(!dLinkVerticesFlag)
		return; //already done prior to here!!!

	const DKinFitParticle* locKinFitParticle;
	deque<const DKinFitParticle*> locConstrainedDecayingParticles;
	if(dDebugLevel > 10)
		cout << "# candidate vertex fits = " << locFinalKinFitParticles_Vertices.size() << endl;
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
							if(dDebugLevel > 10)
								cout << "final constraining: decaying particle constrained by another constraint" << endl;
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
								if(dDebugLevel > 10)
									cout << "init constraining: decaying particle constrained by another constraint" << endl;
								++locNumConstrainedParticles; //decaying particle constrained by another constraint
								break;
							}
						}
					}
				}
			}

			//if enough constrained particles, this fit is OK
			if(dDebugLevel > 10)
				cout << "# vertex constraint constraining particles = " << locNumConstrainedParticles << endl;
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
				if(!locInitialKinFitParticles_Vertices[loc_i].empty())
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
				if(dDebugLevel > 10)
					cout << "vertex fit is OK" << endl;
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

bool DKinFitResults_factory::Calc_VertexGuess(JEventLoop* locEventLoop, const deque<const DKinFitParticle*>& locInitialKinFitParticles, const deque<const DKinFitParticle*>& locFinalKinFitParticles, TVector3& locVertexGuess, deque<const DKinFitParticle*>& locFinalKinFitParticles_Vertex_Updated, int locVertexFindFlag, const map<const DKinFitParticle*, TVector3>& locP4GuessMap, map<const DKinFitParticle*, TVector3>& locDecayingParticleTrackPointGuesses)
{
	//if locVertexFindFlag =: 0 -> 2+ charged tracks else false, 1 -> use decaying info if needed & available else false, 2 -> do best you can with limited info
	if(dDebugLevel > 10)
	{
		cout << "INIT particles in calc vertex guess (mass, q): " << endl;
		for(size_t loc_i = 0; loc_i < locInitialKinFitParticles.size(); ++loc_i)
			cout << locInitialKinFitParticles[loc_i]->Get_Mass() << ", " << locInitialKinFitParticles[loc_i]->Get_Charge() << endl;
		cout << "FINAL particles in calc vertex guess (mass, q): " << endl;
		for(size_t loc_i = 0; loc_i < locFinalKinFitParticles.size(); ++loc_i)
			cout << locFinalKinFitParticles[loc_i]->Get_Mass() << ", " << locFinalKinFitParticles[loc_i]->Get_Charge() << endl;
	}

	//choose midpoint of DOCA line as init guess for vertex
	//do vertex-only kinfit to find the super-init kinfit guess
	deque<const DKinFitParticle*> locVertexFindParticles, locDecayingParticles;
	DKinFitParticleType locKinFitParticleType;
	for(size_t loc_i = 0; loc_i < locFinalKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = locFinalKinFitParticles[loc_i]->Get_KinFitParticleType();
		if((locKinFitParticleType == d_DetectedParticle) && (locFinalKinFitParticles[loc_i]->Get_Charge() != 0))
			locVertexFindParticles.push_back(locFinalKinFitParticles[loc_i]);
		if(locKinFitParticleType == d_DecayingParticle)
			locDecayingParticles.push_back(locFinalKinFitParticles[loc_i]);
	}
	for(size_t loc_i = 0; loc_i < locInitialKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = locInitialKinFitParticles[loc_i]->Get_KinFitParticleType();
		if(locKinFitParticleType == d_DecayingParticle)
			locDecayingParticles.push_back(locInitialKinFitParticles[loc_i]);
	}
	if(dDebugLevel > 50)
	{
		cout << "Vertex-find-particles in calc vertex guess (mass, q): " << endl;
		for(size_t loc_i = 0; loc_i < locVertexFindParticles.size(); ++loc_i)
			cout << locVertexFindParticles[loc_i]->Get_Mass() << ", " << locVertexFindParticles[loc_i]->Get_Charge() << endl;
		cout << "Decaying particles in calc vertex guess (mass, q): " << endl;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
			cout << locDecayingParticles[loc_i]->Get_Mass() << ", " << locDecayingParticles[loc_i]->Get_Charge() << endl;
	}

	deque<const DChargedTrackHypothesis*> locVertexFindParticles_ChargedHypotheses;
	for(size_t loc_i = 0; loc_i < locVertexFindParticles.size(); ++loc_i)
	{
		const DKinematicData* locKinematicData = dKinFitter.Get_Source_FromInput(locVertexFindParticles[loc_i]);
		locVertexFindParticles_ChargedHypotheses.push_back(dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData));
	}

	//if > 2 detected charged particles, don't need to use decaying particles for the guess
	if(locVertexFindParticles.size() >= 2)
	{
		if(dDebugLevel > 20)
			cout << "> 2 charged tracks" << endl;

		//use DReferenceTrajectory::IntersectTracks() to get a very-initial vertex guess, then kinematic fit to get the final vertex guess
		DVector3 locTempInitVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexFindParticles_ChargedHypotheses);
		locVertexGuess.SetXYZ(locTempInitVertex.X(), locTempInitVertex.Y(), locTempInitVertex.Z());
		if(dDebugLevel > 20)
			cout << "init vertex guess from DAnalysisUtilities (Vertex-find-particles) = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;

		//setup fitter
		dKinFitter.Reset_NewFit();
		const DKinFitConstraint_Vertex* locVertexConstraint = dKinFitter.Add_VertexConstraint(locInitialKinFitParticles, locFinalKinFitParticles, locVertexGuess);

		//perform fit
		if(dDebugLevel > 0)
			cout << "Perform init vertex guess kinematic fit" << endl;
		bool locFitStatus = dKinFitter.Fit_Reaction();
		if(dDebugLevel > 0)
			cout << "Init vertex guess kinematic fit finished, status = " << locFitStatus << endl;
		if(locFitStatus)
			locVertexGuess = locVertexConstraint->Get_CommonVertex();
		if(dDebugLevel > 20)
			cout << "final vert guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;

		//save vertex guess as a track point for decaying particles
		map<const DKinFitParticle*, TVector3>::iterator locMapIterator;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
			locDecayingParticleTrackPointGuesses[locDecayingParticles[loc_i]] = locVertexGuess;

		return true;
	}

	//need to use decaying particle info
	deque<const DKinFitParticle*> locVertexFindParticles_Decaying;
	//see how many decaying particles have a known track point (from a previous vertex guess)
	for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	{
		if(locDecayingParticleTrackPointGuesses.find(locDecayingParticles[loc_i]) != locDecayingParticleTrackPointGuesses.end())
		{
			locVertexFindParticles_Decaying.push_back(locDecayingParticles[loc_i]); //point already determined, can use!
			TVector3 locTempDecayingVertexPoint = locDecayingParticleTrackPointGuesses[locDecayingParticles[loc_i]];
			if(dDebugLevel > 15)
				cout << "decaying particle (q, mass = " << locDecayingParticles[loc_i]->Get_Charge() << ", " << locDecayingParticles[loc_i]->Get_Mass() << ") has known point at: " << locTempDecayingVertexPoint.X() << ", " << locTempDecayingVertexPoint.Y() << ", " << locTempDecayingVertexPoint.Z() << endl;
		}
	}
	if(dDebugLevel > 15)
		cout << "#charged, decaying vertex-find-particles = " << locVertexFindParticles.size() << ", " << locVertexFindParticles_Decaying.size() << endl;
	if((locVertexFindParticles.size() + locVertexFindParticles_Decaying.size()) < 2) //still not enough tracks
		return ((locVertexFindFlag == 2) ? true : false);

	//do lazy method: assume no b-field
	deque<DKinematicData*> locVertexFindParticles_Decaying_KinematicData;
	for(size_t loc_i = 0; loc_i < locVertexFindParticles_Decaying.size(); ++loc_i)
	{
		DKinematicData* locDecayingParticleData = new DKinematicData();
		const DKinFitParticle* locDecayingParticle_KinFit = locVertexFindParticles_Decaying[loc_i];
		TVector3 locTempVertex = locDecayingParticleTrackPointGuesses[locDecayingParticle_KinFit];
		locDecayingParticleData->setPosition(DVector3(locTempVertex.X(), locTempVertex.Y(), locTempVertex.Z()));
		TVector3 locTempP3 = locP4GuessMap.find(locDecayingParticle_KinFit)->second;
		locDecayingParticleData->setMomentum(DVector3(locTempP3.X(), locTempP3.Y(), locTempP3.Z()));
		locVertexFindParticles_Decaying_KinematicData.push_back(locDecayingParticleData);
	}
	if(locVertexFindParticles.size() == 1)
	{
		if(dDebugLevel > 20)
			cout << "1 charged track" << endl;

		DVector3 locTempInitVertex = dAnalysisUtilities->Calc_CrudeVertex(locVertexFindParticles_ChargedHypotheses[0], locVertexFindParticles_Decaying_KinematicData);
		for(size_t loc_i = 0; loc_i < locVertexFindParticles_Decaying_KinematicData.size(); ++loc_i)
			delete locVertexFindParticles_Decaying_KinematicData[loc_i]; //delete temp data
		locVertexGuess.SetXYZ(locTempInitVertex.X(), locTempInitVertex.Y(), locTempInitVertex.Z());
		if(dDebugLevel > 20)
			cout << "final vert guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;
		map<const DKinFitParticle*, TVector3>::iterator locMapIterator;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
		{
			if(locDecayingParticleTrackPointGuesses.find(locDecayingParticles[loc_i]) == locDecayingParticleTrackPointGuesses.end())
				locDecayingParticleTrackPointGuesses[locDecayingParticles[loc_i]] = locVertexGuess; //save point for decaying tracks (if not done so already
		}
		return true;
	}
	//FINISH ME: TWO DECAYING PARTICLES
	return true;
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

void DKinFitResults_factory::Build_KinFitResults(const DParticleCombo* locParticleCombo, const map<const DKinFitParticle*, pair<Particle_t, deque<const DKinematicData*> > >& locInitDecayingKinFitParticles)
{
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	DKinFitResults* locKinFitResults = new DKinFitResults();
	locKinFitResults->Set_ParticleCombo(locParticleCombo);

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

bool DKinFitResults_factory::Find_ConstrainableParticles(deque<pair<const DKinFitParticle*, size_t> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_P4s, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_P4s, const deque<size_t>& locConstraintsSetIndices)
{
	DKinFitParticleType locKinFitParticleType;
	for(size_t loc_i = 0; loc_i < locInitialKinFitParticles_P4s.size(); ++loc_i)
	{
		bool locConstraintAlreadySetFlag = false;
		for(size_t loc_j = 0; loc_j < locConstraintsSetIndices.size(); ++loc_j)
		{
			if(locConstraintsSetIndices[loc_j] != loc_i)
				continue;
			locConstraintAlreadySetFlag = true;
			break;
		}
		if(locConstraintAlreadySetFlag)
			continue; //constraint already set
		size_t locNumParticlesNeedToBeConstrained = 0;
		const DKinFitParticle* locConstrainableParticle = NULL;
		for(size_t loc_j = 0; loc_j < locInitialKinFitParticles_P4s[loc_i].size(); ++loc_j)
		{
			locKinFitParticleType = locInitialKinFitParticles_P4s[loc_i][loc_j]->Get_KinFitParticleType();
			if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
				continue;
			bool locParticleAlreadyConstrainedFlag = false;
			for(size_t loc_k = 0; loc_k < locConstrainedParticles.size(); ++loc_k)
			{
				if(locConstrainedParticles[loc_k] != locInitialKinFitParticles_P4s[loc_i][loc_j])
					continue;
				locParticleAlreadyConstrainedFlag = true;
				break;
			}
			if(locParticleAlreadyConstrainedFlag)
				continue;
			++locNumParticlesNeedToBeConstrained;
			locConstrainableParticle = locInitialKinFitParticles_P4s[loc_i][loc_j];
		}
		for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_P4s[loc_i].size(); ++loc_j)
		{
			locKinFitParticleType = locFinalKinFitParticles_P4s[loc_i][loc_j]->Get_KinFitParticleType();
			if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
				continue;
			bool locParticleAlreadyConstrainedFlag = false;
			for(size_t loc_k = 0; loc_k < locConstrainedParticles.size(); ++loc_k)
			{
				if(locConstrainedParticles[loc_k] != locFinalKinFitParticles_P4s[loc_i][loc_j])
					continue;
				locParticleAlreadyConstrainedFlag = true;
				break;
			}
			if(locParticleAlreadyConstrainedFlag)
				continue;
			++locNumParticlesNeedToBeConstrained;
			locConstrainableParticle = locFinalKinFitParticles_P4s[loc_i][loc_j];
		}
		if(locNumParticlesNeedToBeConstrained == 1) //else too many unconstrained particles in it's step to be able to constrain it right away!!
			locConstrainableParticles.push_back(pair<const DKinFitParticle*, size_t>(locConstrainableParticle, loc_i));
	}
	return (!locConstrainableParticles.empty());
}

void DKinFitResults_factory::Calc_P4Guess(pair<const DKinFitParticle*, size_t>& locConstrainableParticle, const deque<deque<const DKinFitParticle*> >& locInitialKinFitParticles_P4s, const deque<deque<const DKinFitParticle*> >& locFinalKinFitParticles_P4s, TVector3& locMomentum)
{
	const DKinFitParticle* locParticleToConstrain = locConstrainableParticle.first;
	size_t locConstraintIndex = locConstrainableParticle.second;
	const DKinFitParticle* locKinFitParticle;

	if(dDebugLevel > 5)
		cout << "particle to constrain q, mass = " << locParticleToConstrain->Get_Charge() << ", " << locParticleToConstrain->Get_Mass() << endl;

	bool locConstrainedParticleIsInInitialState = false;
	for(size_t loc_j = 0; loc_j < locInitialKinFitParticles_P4s[locConstraintIndex].size(); ++loc_j)
	{
		locKinFitParticle = locInitialKinFitParticles_P4s[locConstraintIndex][loc_j];
		if(locKinFitParticle == locParticleToConstrain)
			locConstrainedParticleIsInInitialState = true;
		else
		{
			TVector3 locParticleMomentum = locKinFitParticle->Get_Momentum();
			if(dDebugLevel > 15)
				cout << "initial particle mass, q, pxyz = " << locKinFitParticle->Get_Mass() << ", " << locKinFitParticle->Get_Charge() << ", " << locParticleMomentum.Px() << ", " << locParticleMomentum.Py() << ", " << locParticleMomentum.Pz() << endl;
			locMomentum += locParticleMomentum;
		}
	}
	for(size_t loc_j = 0; loc_j < locFinalKinFitParticles_P4s[locConstraintIndex].size(); ++loc_j)
	{
		locKinFitParticle = locFinalKinFitParticles_P4s[locConstraintIndex][loc_j];
		if(locKinFitParticle != locParticleToConstrain)
		{
			TVector3 locParticleMomentum = locKinFitParticle->Get_Momentum();
			locMomentum -= locParticleMomentum;
			if(dDebugLevel > 15)
				cout << "final particle mass, q, pxyz = " << locKinFitParticle->Get_Mass() << ", " << locKinFitParticle->Get_Charge() << ", " << locParticleMomentum.Px() << ", " << locParticleMomentum.Py() << ", " << locParticleMomentum.Pz() << endl;
		}
	}

	if(locConstrainedParticleIsInInitialState)
		locMomentum *= -1.0;

	if(dDebugLevel > 15)
		cout << "p3 guess = " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << endl;
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

