#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DKinFitUtils_GlueX.h"

/******************************************************************** INITIALIZE *******************************************************************/

DKinFitUtils_GlueX::DKinFitUtils_GlueX(const DMagneticFieldMap* locMagneticFieldMap, const DAnalysisUtilities* locAnalysisUtilities) : 
dMagneticFieldMap(locMagneticFieldMap), dAnalysisUtilities(locAnalysisUtilities)
{
	dIncludeBeamlineInVertexFitFlag = false;
	dWillBeamHaveErrorsFlag = false; //Until fixed!
	dEventNumber = 0;

	//fill buffer: the larger the number, the more memory it takes. the smaller, the more locking is needed
	dNumFillBufferMatrices = 10;
	dNumFillBufferParticles = 10;
	dTargetMaxNumAvailableParticles = 250000;

	dApplication = dynamic_cast<DApplication*>(japp);
	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);
}

DKinFitUtils_GlueX::DKinFitUtils_GlueX(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	dApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	dMagneticFieldMap = dApplication->GetBfield(locEventLoop->GetJEvent().GetRunNumber());

	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);
	dWillBeamHaveErrorsFlag = false; //Until fixed!
	dIncludeBeamlineInVertexFitFlag = false;

	dEventNumber = locEventLoop->GetJEvent().GetEventNumber();

	//fill buffer: the larger the number, the more memory it takes. the smaller, the more locking is needed
	dNumFillBufferMatrices = 10;
	dNumFillBufferParticles = 10;
	dTargetMaxNumAvailableParticles = 250000;
}

void DKinFitUtils_GlueX::Set_MaxPoolSizes(size_t locNumReactions, size_t locExpectedNumCombos)
{
	//final x2: input/output
	Set_MaxKinFitParticlePoolSize(10*locNumReactions*locExpectedNumCombos*2);

	Set_MaxKinFitConstraintVertexPoolSize(2*2*locNumReactions*locExpectedNumCombos*2); //extra x2: guess fit!
	Set_MaxKinFitConstraintSpacetimePoolSize(2*2*locNumReactions*locExpectedNumCombos*2); //extra x2: guess fit!
	Set_MaxKinFitConstraintP4PoolSize(locNumReactions*locExpectedNumCombos*2);
	Set_MaxKinFitConstraintMassPoolSize(2*locNumReactions*locExpectedNumCombos*2);

	Set_MaxKinFitChainPoolSize(locNumReactions*locExpectedNumCombos);
	Set_MaxKinFitChainStepPoolSize(3*locNumReactions*locExpectedNumCombos);

	Set_MaxSymMatrixPoolSize(10*locNumReactions*locExpectedNumCombos*2);
	dTargetMaxNumAvailableParticles = 5000*locNumReactions;
}

/*********************************************************** OVERRIDE BASE CLASS FUNCTIONS *********************************************************/

void DKinFitUtils_GlueX::Reset_NewEvent(uint64_t locEventNumber)
{
	dEventNumber = locEventNumber;
	Reset_NewEvent();
}

void DKinFitUtils_GlueX::Reset_NewEvent(void)
{
	dParticleMap_SourceToInput_Beam.clear();
	dParticleMap_SourceToInput_DetectedParticle.clear();
	dParticleMap_SourceToInput_Shower.clear();
	dParticleMap_SourceToInput_Target.clear();
	dParticleMap_SourceToInput_Decaying.clear();
	dParticleMap_SourceToInput_Missing.clear();

	dParticleMap_InputToSource_JObject.clear();
	dParticleMap_InputToSource_Decaying.clear();

	Reset_ParticleMemory();
	DKinFitUtils::Reset_NewEvent();
}

bool DKinFitUtils_GlueX::Get_IncludeBeamlineInVertexFitFlag(void) const
{
	return dIncludeBeamlineInVertexFitFlag; //at least until covariance matrix is set for beam photons
}

bool DKinFitUtils_GlueX::Get_IsDetachedVertex(int locPDG_PID) const
{
	return IsDetachedVertex(PDGtoPType(locPDG_PID));
}

bool DKinFitUtils_GlueX::Get_IsBFieldNearBeamline(void) const
{
	if(dMagneticFieldMap == NULL)
		return false;

	return (dynamic_cast<const DMagneticFieldMapNoField*>(dMagneticFieldMap) == NULL);
}

TVector3 DKinFitUtils_GlueX::Get_BField(const TVector3& locPosition) const
{
	if(dMagneticFieldMap == NULL)
		return TVector3(0.0, 0.0, 0.0);

	double locBx, locBy, locBz;
	dMagneticFieldMap->GetField(locPosition.X(), locPosition.Y(), locPosition.Z(), locBx, locBy, locBz);
	return (TVector3(locBx, locBy, locBz));
}

/****************************************************************** MANAGE MEMORY ******************************************************************/

deque<DKinFitParticle*>& DKinFitUtils_GlueX::Get_AvailableParticleDeque(void) const
{
	//static: shared amongst all threads
	//Must call this function within a lock!!
	static deque<DKinFitParticle*> locAvailableParticles;
	return locAvailableParticles;
}

DKinFitParticle* DKinFitUtils_GlueX::Get_KinFitParticleResource(void)
{
	//if kinfit pool (buffer) is empty, use shared pool to retrieve a new batch of particles
	if(Get_KinFitParticlePoolAvailableSize() == 0)
		Acquire_Particles(dNumFillBufferParticles);

	return DKinFitUtils::Get_KinFitParticleResource();
}

void DKinFitUtils_GlueX::Reset_ParticleMemory(void)
{
	//Access combo resource pool
	bool locDeleteParticlesFlag = false;
	japp->WriteLock("DKinFitParticle_Memory"); //LOCK
	{
		deque<DKinFitParticle*>& locAvailableParticles = Get_AvailableParticleDeque();

		//Memory fragmentation seems to be a very big problem, and these objects use the most memory
		//So, don't delete them. To delete them, just uncomment the below lines.

		//if available > max, wipe all used
//		locDeleteParticlesFlag = (locAvailableParticles.size() >= dTargetMaxNumAvailableParticles);
//		if(!locDeleteParticlesFlag)
			std::move(dKinFitParticlePool_Acquired.begin(), dKinFitParticlePool_Acquired.end(), std::back_inserter(locAvailableParticles));
	}
	japp->Unlock("DKinFitParticle_Memory"); //UNLOCK

	//delete combos if necessary
	if(locDeleteParticlesFlag)
	{
		for(auto& locParticle : dKinFitParticlePool_Acquired)
			delete locParticle;
	}

	//clear thread-local pool
	dKinFitParticlePool_Acquired.clear();
}

void DKinFitUtils_GlueX::Acquire_Particles(size_t locNumRequestedParticles)
{
	//We must have the correct event number, so that we know when it's safe to recycle the memory for the next event.
	deque<DKinFitParticle*> locAcquiredParticles;

	//Access resource pool
	japp->WriteLock("DKinFitParticle_Memory"); //LOCK
	{
		deque<DKinFitParticle*>& locAvailableParticles = Get_AvailableParticleDeque();

		//Get resources if available
		if(locAvailableParticles.size() <= locNumRequestedParticles)
		{
			//Move the whole deque
			std::move(locAvailableParticles.begin(), locAvailableParticles.end(), std::back_inserter(locAcquiredParticles));
			locAvailableParticles.clear();
		}
		else //Move the desired range
		{
			size_t locNewSize = locAvailableParticles.size() - locNumRequestedParticles;
			auto locMoveEdgeIterator = std::next(locAvailableParticles.rbegin(), locNumRequestedParticles);
			std::move(locAvailableParticles.rbegin(), locMoveEdgeIterator, std::back_inserter(locAcquiredParticles));
			locAvailableParticles.resize(locNewSize);
		}
	}
	japp->Unlock("DKinFitParticle_Memory"); //UNLOCK

	//set matrix pointer as null (necessary before "recycling" below)
		//Recycle_Particles() also recycles matrix memory, want to avoid that, so null them first
	for(auto& locParticle : locAcquiredParticles)
		locParticle->Set_CovarianceMatrix(nullptr);

	//make new particles if necessary
	while(locAcquiredParticles.size() < locNumRequestedParticles)
		locAcquiredParticles.push_back(new DKinFitParticle());

	//Store the acquired particles to the dKinFitParticlePool_Available buffer by "recycling" them
		//these only live in the "available" pool, and aren't set in the "all" pool
		//when the pools are reset for a new event, the buffer is cleared and the utils forget all about them
		//thus, the memory is managed by DKinFitUtils_GlueX, and not by DKinFitUtils
	set<DKinFitParticle*> locParticlesToRecycle(locAcquiredParticles.begin(), locAcquiredParticles.end());
	Recycle_Particles(locParticlesToRecycle);

	//Register as acquired-by this thread
	std::move(locAcquiredParticles.begin(), locAcquiredParticles.end(), std::back_inserter(dKinFitParticlePool_Acquired));
}

size_t DKinFitUtils_GlueX::Get_KinFitParticlePoolSize_Shared(void) const
{
	size_t locNumParticles = 0;
	
	//Access resource pool
	japp->WriteLock("DKinFitParticle_Memory"); //LOCK
	{
		locNumParticles = Get_AvailableParticleDeque().size();
	}
	japp->Unlock("DKinFitParticle_Memory"); //UNLOCK

	return locNumParticles;
}

void DKinFitUtils_GlueX::Recycle_DetectedDecayingParticles(map<DKinFitParticle*, DKinFitParticle*>& locDecayingToDetectedParticleMap)
{
	set<DKinFitParticle*> locParticlesToRecycle;
	for(auto& locParticlePair : locDecayingToDetectedParticleMap)
		locParticlesToRecycle.insert(locParticlePair.second);

	Recycle_Particles(locParticlesToRecycle);
	locDecayingToDetectedParticleMap.clear();
}

/****************************************************************** MAKE PARTICLES *****************************************************************/

DKinFitParticle* DKinFitUtils_GlueX::Make_BeamParticle(const DBeamPhoton* locBeamPhoton)
{
	pair<const DBeamPhoton*, const DEventRFBunch*> locSourcePair(locBeamPhoton, NULL);
	if(dParticleMap_SourceToInput_Beam.find(locSourcePair) != dParticleMap_SourceToInput_Beam.end())
		return dParticleMap_SourceToInput_Beam[locSourcePair]; //not unique, return existing

	TLorentzVector locSpacetimeVertex(Make_TVector3(locBeamPhoton->position()), locBeamPhoton->time());
	TVector3 locMomentum = Make_TVector3(locBeamPhoton->momentum());
	Particle_t locPID = locBeamPhoton->PID();

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_BeamParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), 
		locSpacetimeVertex, locMomentum, locBeamPhoton->errorMatrix());
	dParticleMap_SourceToInput_Beam[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locBeamPhoton;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_BeamParticle(const DBeamPhoton* locBeamPhoton, const DEventRFBunch* locEventRFBunch)
{
	pair<const DBeamPhoton*, const DEventRFBunch*> locSourcePair(locBeamPhoton, locEventRFBunch);
	if(dParticleMap_SourceToInput_Beam.find(locSourcePair) != dParticleMap_SourceToInput_Beam.end())
		return dParticleMap_SourceToInput_Beam[locSourcePair]; //not unique, return existing

	//set rf time for beam particle
	TLorentzVector locSpacetimeVertex(Make_TVector3(locBeamPhoton->position()), locEventRFBunch->dTime);
	TVector3 locMomentum = Make_TVector3(locBeamPhoton->momentum());
	Particle_t locPID = locBeamPhoton->PID();

	//set rf time variance in covariance matrix
	TMatrixFSym* locCovarianceMatrix = Get_SymMatrixResource(7);
	*locCovarianceMatrix = *(locBeamPhoton->errorMatrix());
	(*locCovarianceMatrix)(6, 6) = locEventRFBunch->dTimeVariance;
	//zero the correlation terms
	for(int loc_i = 0; loc_i < 6; ++loc_i)
	{
		(*locCovarianceMatrix)(6, loc_i) = 0.0;
		(*locCovarianceMatrix)(loc_i, 6) = 0.0;
	}

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_BeamParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), 
		locSpacetimeVertex, locMomentum, locCovarianceMatrix);
	dParticleMap_SourceToInput_Beam[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locBeamPhoton;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_DetectedParticle(const DKinematicData* locKinematicData)
{
	if(dParticleMap_SourceToInput_DetectedParticle.find(locKinematicData) != dParticleMap_SourceToInput_DetectedParticle.end())
		return dParticleMap_SourceToInput_DetectedParticle[locKinematicData]; //not unique, return existing

	TLorentzVector locSpacetimeVertex(Make_TVector3(locKinematicData->position()), locKinematicData->time());
	TVector3 locMomentum = Make_TVector3(locKinematicData->momentum());
	Particle_t locPID = locKinematicData->PID();

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_DetectedParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), 
		locSpacetimeVertex, locMomentum, locKinematicData->errorMatrix());
	dParticleMap_SourceToInput_DetectedParticle[locKinematicData] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locKinematicData;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_DetectedParticle(DKinFitParticle* locDecayingKinFitParticle)
{
	if(locDecayingKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
	{
		cout << "WARNING: NON-DECAYING PARTICLE AS INPUT TO DKinFitUtils_GlueX::Make_DetectedParticle(). RETURNING NULL." << endl;
		return NULL;
	}

	DKinFitParticle* locDetectedKinFitParticle = DKinFitUtils::Make_DetectedParticle(locDecayingKinFitParticle->Get_PID(), 
		locDecayingKinFitParticle->Get_Charge(), locDecayingKinFitParticle->Get_Mass(), locDecayingKinFitParticle->Get_SpacetimeVertex(), 
		locDecayingKinFitParticle->Get_Momentum(), Clone_SymMatrix(locDecayingKinFitParticle->Get_CovarianceMatrix()));

	return locDetectedKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_DetectedShower(const DNeutralShower* locNeutralShower, Particle_t locPID)
{
	pair<const DNeutralShower*, Particle_t> locSourcePair(locNeutralShower, locPID);
	if(dParticleMap_SourceToInput_Shower.find(locSourcePair) != dParticleMap_SourceToInput_Shower.end())
		return dParticleMap_SourceToInput_Shower[locSourcePair]; //not unique, return existing

	//use DNeutralShower object (doesn't make assumption about vertex!)
	TLorentzVector locShowerSpacetime = Make_TLorentzVector(locNeutralShower->dSpacetimeVertex);
	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_DetectedShower(PDGtype(locPID), ParticleMass(locPID), locShowerSpacetime, 
		locNeutralShower->dEnergy, &locNeutralShower->dCovarianceMatrix);

	dParticleMap_SourceToInput_Shower[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locNeutralShower;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_TargetParticle(Particle_t locPID)
{
	if(dParticleMap_SourceToInput_Target.find(locPID) != dParticleMap_SourceToInput_Target.end())
		return dParticleMap_SourceToInput_Target[locPID]; //not unique, return existing

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_TargetParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMap_SourceToInput_Target[locPID] = locKinFitParticle;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_MissingParticle(Particle_t locPID)
{
	if(dParticleMap_SourceToInput_Missing.find(locPID) != dParticleMap_SourceToInput_Missing.end())
		return dParticleMap_SourceToInput_Missing[locPID]; //not unique, return existing

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_MissingParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMap_SourceToInput_Missing[locPID] = locKinFitParticle;
	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils_GlueX::Make_DecayingParticle(Particle_t locPID, const set<DKinFitParticle*>& locFromInitialState, const set<DKinFitParticle*>& locFromFinalState)
{
	DDecayingParticleInfo locDecayingParticleInfo(locPID, locFromInitialState, locFromFinalState);
	if(dParticleMap_SourceToInput_Decaying.find(locDecayingParticleInfo) != dParticleMap_SourceToInput_Decaying.end())
		return dParticleMap_SourceToInput_Decaying[locDecayingParticleInfo]; //not unique, return existing

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_DecayingParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), locFromInitialState, locFromFinalState);
	dParticleMap_SourceToInput_Decaying[locDecayingParticleInfo] = locKinFitParticle;
	dParticleMap_InputToSource_Decaying.insert(pair<DKinFitParticle*, DDecayingParticleInfo>(locKinFitParticle, locDecayingParticleInfo));
	return locKinFitParticle;
}

/**************************************************************** MAKE DKINFITCHAIN ****************************************************************/

const DKinFitChain* DKinFitUtils_GlueX::Make_KinFitChain(const DParticleCombo* locParticleCombo, DKinFitType locKinFitType)
{
	//locKinFitType input in case want to do a different fit
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: Create DKinFitChain." << endl;

	DKinFitChain* locKinFitChain = Get_KinFitChainResource();
	locKinFitChain->Set_DefinedParticleStepIndex(-1); //unless changed below

	//Make chain, excluding decaying particles
		//They must be created using the detected particles, so just create those first
	map<size_t, size_t> locStepCreationMap; //key is source particle combo step, value is created kinfit step index
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		if(locStepCreationMap.find(loc_i) != locStepCreationMap.end())
			continue; //already did this step

		//Start and register new step
		DKinFitChainStep* locKinFitChainStep = Get_KinFitChainStepResource();
		locKinFitChainStep->Set_ConstrainDecayingMassFlag(false); //unless changed later

		//Create it
		Make_KinFitChainStep(locParticleCombo, locKinFitType, loc_i, locKinFitChain, locKinFitChainStep, locStepCreationMap);
		locKinFitChain->Add_KinFitChainStep(locKinFitChainStep);
	}

	//Now define the decaying particles using the detected & beam particles

	//start looping from the back, using invariant mass wherever possible
		//not possible if missing decay product: will then compute via missing mass in the next step
	//but first, figure out which steps we must use missing mass for (and thus must skip on this pass below)
	set<size_t> locMissingMassSteps;
	int locCurrentStepIndex = locKinFitChain->Get_DefinedParticleStepIndex();
	while(locCurrentStepIndex != -1)
	{
		if((locCurrentStepIndex == 0) && locKinFitChain->Get_KinFitChainStep(0)->Get_InitialParticles().empty())
			break; //is an open-ended decaying particle in the initial state: is ok to define via invariant mass
		locMissingMassSteps.insert(locCurrentStepIndex);
		locCurrentStepIndex = locKinFitChain->Get_KinFitChainStep(locCurrentStepIndex)->Get_InitialParticleDecayFromStepIndex();
	}

	//now do the invariant mass loop
	for(int loc_i = locParticleCombo->Get_NumParticleComboSteps() - 1; loc_i >= 0; --loc_i)
	{
		//get steps
		int locKinFitStepIndex = locStepCreationMap[loc_i];
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		DKinFitChainStep* locKinFitChainStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locKinFitStepIndex));

		//skip steps that must defined via missing mass
		if(locMissingMassSteps.find(locKinFitStepIndex) != locMissingMassSteps.end())
			continue; //decay products contain a missing or open-ended-decaying particle

		//check if initial particle already created for this step (i.e. beam)
		if(!locKinFitChainStep->Get_InitialParticles().empty())
			continue; //only do decaying particles

		//don't create particle for omega, etc. (will confuse kinfitter)
		Particle_t locPID = locParticleComboStep->Get_InitialParticleID();
		if(!IsFixedMass(locPID))
			continue;

		//since going in reverse order, all decay products are ready: create the decaying particle
		set<DKinFitParticle*> locFromInitialState; //empty
		DKinFitParticle* locDecayingParticle = Make_DecayingParticle(locPID, locFromInitialState, locKinFitChainStep->Get_FinalParticles());

		//set decaying particle in the chain
		locKinFitChainStep->Add_InitialParticle(locDecayingParticle);
		int locProducedStepIndex = locStepCreationMap[locParticleComboStep->Get_InitialParticleDecayFromStepIndex()];
		if((locProducedStepIndex < 0) || (locProducedStepIndex == loc_i))
			continue; //intial-state open-ended decaying particle

		DKinFitChainStep* locProductionStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locProducedStepIndex));
		locProductionStep->Add_FinalParticle(locDecayingParticle);
		locKinFitChain->Set_DecayStepIndex(locDecayingParticle, locKinFitStepIndex);

		//mark the progress, and go to the next step
		locStepCreationMap.erase(loc_i);
	}

	//now loop from the front, using missing mass for the rest
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locStepCreationMap.find(loc_i) == locStepCreationMap.end())
			continue;
		int locKinFitStepIndex = locStepCreationMap[loc_i];
		DKinFitChainStep* locKinFitChainStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locKinFitStepIndex));

		//create decaying particle in final state: first figure out which one needs a particle
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			//DecayStepIndex: >= 0 if decaying, where the # is the step representing the particle decay
			int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
			if(locDecayStepIndex < 0)
				continue; //not decaying

			int locKinFitDecayStepIndex = locStepCreationMap[locDecayStepIndex];
			if(!locKinFitChain->Get_KinFitChainStep(locKinFitDecayStepIndex)->Get_InitialParticles().empty())
				continue; //decaying particle already created

			//don't create particle for omega, etc. (will confuse kinfitter)
			Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);
			if(!IsFixedMass(locPID))
				continue;

			//make decaying particle
			DKinFitParticle* locDecayingParticle = Make_DecayingParticle(locPID, locKinFitChainStep->Get_InitialParticles(), locKinFitChainStep->Get_FinalParticles());
			locKinFitChainStep->Add_FinalParticle(locDecayingParticle);

			DKinFitChainStep* locDecayStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locKinFitDecayStepIndex));
			locDecayStep->Add_InitialParticle(locDecayingParticle);
			locKinFitChain->Set_DecayStepIndex(locDecayingParticle, locKinFitDecayStepIndex);
		}
	}

	if(dDebugLevel > 10)
	{
		cout << "DKinFitUtils_GlueX: DKinFitChain Created. Printing:" << endl;
		locKinFitChain->Print_InfoToScreen();
	}

	return locKinFitChain;
}

void DKinFitUtils_GlueX::Make_KinFitChainStep(const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, size_t locStepIndex, DKinFitChain* locKinFitChain, DKinFitChainStep* locKinFitChainStep, map<size_t, size_t>& locStepCreationMap)
{
	//get the steps
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	const DReaction* locReaction = locParticleCombo->Get_Reaction();
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	int locKinFitStepIndex = locKinFitChain->Get_NumKinFitChainSteps();

	//if doing a vertex fit, see which neutral particles can be treated as showers
	set<pair<int, int> > locKinFitVertexParticles;
	bool locSpactimeIsFitFlag = (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit);
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
		locKinFitVertexParticles = Get_KinFitVertexParticles(locReaction);

	//initial particle
	if(locKinFitChainStep->Get_InitialParticles().empty()) //else step already started: initial particle doesn't have fixed mass
	{
		locKinFitChainStep->Set_InitialParticleDecayFromStepIndex(-1); //unless set otherwise below (enclosed decaying particle)
		//initial beam particle
		const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locParticleComboStep->Get_InitialParticle());
		if(locBeamPhoton != NULL)
			locKinFitChainStep->Add_InitialParticle(Make_BeamParticle(locBeamPhoton));
		else //decaying particle
		{
			if(locStepIndex == 0) //open-ended
				locKinFitChain->Set_DefinedParticleStepIndex(locKinFitStepIndex);
			else //enclosed
			{
				int locDecayFromStepIndex = locStepCreationMap[locParticleComboStep->Get_InitialParticleDecayFromStepIndex()];
				locKinFitChainStep->Set_InitialParticleDecayFromStepIndex(locDecayFromStepIndex);
				Particle_t locPID = locParticleComboStep->Get_InitialParticleID();
				if(IsFixedMass(locPID))
					locKinFitChainStep->Set_ConstrainDecayingMassFlag(locReactionStep->Get_KinFitConstrainInitMassFlag());
			}
		}
		//skip creating decaying particles for now: will create later

		//target particle
		Particle_t locTargetPID = locParticleComboStep->Get_TargetParticleID();
		if(locTargetPID != Unknown)
			locKinFitChainStep->Add_InitialParticle(Make_TargetParticle(locTargetPID));
	}

	//final state particles
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_j);
		const DKinematicData* locKinematicData = locParticleComboStep->Get_FinalParticle(loc_j);
		Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_j);

		if(locDecayStepIndex == -1) //missing particle
		{
			locKinFitChain->Set_DefinedParticleStepIndex(locKinFitStepIndex);
			if(locPID == Unknown)
			{
				locKinFitChain->Set_IsInclusiveChannelFlag(true);
				continue;
			}
			locKinFitChainStep->Add_FinalParticle(Make_MissingParticle(locPID));
		}
		else if(locDecayStepIndex >= 0) //decaying particle
		{
			if(IsFixedMass(locPID))
				continue; //skip for now, will create later
			//e.g. omega: cannot constrain: add its decay step to this one
			Make_KinFitChainStep(locParticleCombo, locKinFitType, locDecayStepIndex, locKinFitChain, locKinFitChainStep, locStepCreationMap);
		}
		else if(ParticleCharge(locPID) == 0) //detected neutral
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);

			//Determine whether we should use the particle or the shower object
			pair<int, int> locParticlePair(locStepIndex, loc_j);
			bool locNeutralShowerFlag = (locKinFitVertexParticles.find(locParticlePair) != locKinFitVertexParticles.end());
			if((ParticleMass(locPID) > 0.0) && !locSpactimeIsFitFlag)
				locNeutralShowerFlag = false; //massive shower momentum is defined by t, which isn't fit: use particle

			if(!locNeutralShowerFlag)
				locKinFitChainStep->Add_FinalParticle(Make_DetectedParticle(locNeutralParticleHypothesis));
			else //in a vertex constraint: make shower
			{
				const DNeutralShower* locNeutralShower = NULL;
				locNeutralParticleHypothesis->GetSingle(locNeutralShower);
				locKinFitChainStep->Add_FinalParticle(Make_DetectedShower(locNeutralShower, locNeutralParticleHypothesis->PID()));
			}
		}
		else //detected charged track
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locKinematicData);
			locKinFitChainStep->Add_FinalParticle(Make_DetectedParticle(locChargedTrackHypothesis));
		}
	}

	locStepCreationMap.insert(pair<size_t, size_t>(locStepIndex, locKinFitStepIndex));
}

/**************************************************************** MAKE CONSTRAINTS *****************************************************************/

set<DKinFitConstraint*> DKinFitUtils_GlueX::Create_Constraints(const DReactionVertexInfo* locReactionVertexInfo, const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain, DKinFitType locKinFitType, deque<DKinFitConstraint_Vertex*>& locSortedVertexConstraints)
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: Create constraints." << endl;

	//All constraints
	set<DKinFitConstraint*> locAllConstraints;

	//Create Mass Constraints
	set<DKinFitConstraint_Mass*> locMassConstraints;
	map<DKinFitParticle*, DKinFitConstraint_Mass*> locParticleMassConstraintMap;
	map<DKinFitParticle*, size_t> locParticleDecayStepMap; //key is decaying particle, value is step index
	map<size_t, DKinFitConstraint_Mass*> locStepMassConstraintMap;
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//loop over steps, but skip init step (if open-ended decaying, do p4 constraint instead)
		for(size_t loc_i = 1; loc_i < locKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
		{
			const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(loc_i);
			if(!locKinFitChainStep->Get_ConstrainDecayingMassFlag())
				continue; //don't apply mass constraint to this step

			set<DKinFitParticle*> locInitialParticles = locKinFitChainStep->Get_InitialParticles();
			set<DKinFitParticle*>::iterator locParticleIterator = locInitialParticles.begin();
			for(; locParticleIterator != locInitialParticles.end(); ++locParticleIterator)
			{
				if((*locParticleIterator)->Get_KinFitParticleType() != d_DecayingParticle)
					continue; //not a decaying particle
				DKinFitConstraint_Mass* locMassConstraint = Make_MassConstraint(*locParticleIterator);
				locMassConstraints.insert(locMassConstraint);
				locParticleMassConstraintMap[*locParticleIterator] = locMassConstraint;
				locParticleDecayStepMap[*locParticleIterator] = loc_i;
				locStepMassConstraintMap[loc_i] = locMassConstraint;
			}
		}
	}

	//Create P4 Constraint
		//don't do if inclusive reaction
		//pick the step containing the defined (missing or open-ended-decaying) particle
		//if no defined particle, use the first step
	DKinFitConstraint_P4* locP4Constraint = NULL;
	if(!locKinFitChain->Get_IsInclusiveChannelFlag() && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		int locDefinedParticleStepIndex = locKinFitChain->Get_DefinedParticleStepIndex();
		int locP4StepIndex = (locDefinedParticleStepIndex >= 0) ? locDefinedParticleStepIndex : 0;
		const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(locP4StepIndex);
		locP4Constraint = Make_P4Constraint(locKinFitChainStep->Get_InitialParticles(), locKinFitChainStep->Get_FinalParticles());

		//OK, now, check to see if the system is overly constrained: 
			//there must be at least one particle with non-zero errors in the p4 constraint that is NOT in a mass constraint
		bool locNonZeroErrorFlag = false;
		set<size_t> locP4ConstrainedParticleSteps;

		set<DKinFitParticle*> locAllParticles = locP4Constraint->Get_AllParticles();
		set<DKinFitParticle*>::iterator locParticleIterator = locAllParticles.begin();
		for(; locParticleIterator != locAllParticles.end(); ++locParticleIterator)
		{
			DKinFitParticle* locKinFitParticle = *locParticleIterator;
			DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

			//check if decaying particle mass is not constrained
			if(locKinFitParticleType == d_DecayingParticle)
			{
				if(locParticleMassConstraintMap.find(locKinFitParticle) == locParticleMassConstraintMap.end())
				{
					locNonZeroErrorFlag = true; //not constrained: we are OK
					break;
				}
				locP4ConstrainedParticleSteps.insert(locParticleDecayStepMap[locKinFitParticle]);
			}
			if(locKinFitParticle->Get_CovarianceMatrix() == NULL)
				continue;
			if((locKinFitParticleType == d_BeamParticle) && !dWillBeamHaveErrorsFlag)
				continue;

			locNonZeroErrorFlag = true; //this particle has non-zero errors: we are OK
			break;
		}

		if(!locNonZeroErrorFlag)
		{
			//system is over-constrained: we must delete a constraint
			//if there is a missing/open-ended particle: delete a mass constraint; else, delete the p4 constraint
			if((locP4Constraint->Get_DefinedParticle() == NULL) || locP4ConstrainedParticleSteps.empty())
				locP4Constraint = NULL; //remove the p4 constraint
			else //remove a mass constraint: delete the one in the earliest step (so consistent) (will be by missing mass if present)
			{
				size_t locEarliestStepIndex = *locP4ConstrainedParticleSteps.begin();
				locMassConstraints.erase(locStepMassConstraintMap[locEarliestStepIndex]); //remove constraint
			}
		}

		//set init p3 guess
		if(locP4Constraint != NULL)
		{
			if(locP4Constraint->Get_DefinedParticle() != NULL)
			{
				DLorentzVector locDefinedP4;
				if(locP4Constraint->Get_IsDefinedParticleInFinalState())
					locDefinedP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo, false);
				else
					locDefinedP4 = dAnalysisUtilities->Calc_FinalStateP4(locParticleCombo, 0, false);
				TVector3 locInitP3Guess(locDefinedP4.Px(), locDefinedP4.Py(), locDefinedP4.Pz());
				locP4Constraint->Set_InitP3Guess(locInitP3Guess);
			}
			else
				locP4Constraint->Set_InitP3Guess(TVector3(0.0, 0.0, 0.0));
		}
	}

	//Set P4 & Mass constraints
	if(locP4Constraint != NULL)
		locAllConstraints.insert(locP4Constraint);
	locAllConstraints.insert(locMassConstraints.begin(), locMassConstraints.end());

	//Create Vertex Constraints
	//VERTEX OR SPACETIME: Group particles by detached vertex (one deque for each constraint/vertex)
	locSortedVertexConstraints.clear();
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		bool locSpacetimeFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));
		for(locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		{
			auto locDX4 = locParticleCombo->Get_ParticleComboStep(locStepVertexInfo->Get_StepIndices().front())->Get_SpacetimeVertex();
			TLorentzVector locX4(locDX4.X(), locDX4.Y(), locDX4.Z(), locDX4.T());
			if(locSpacetimeFitFlag)
				locSortedVertexConstraints.push_back(Make_SpacetimeConstraint(locStepVertexInfo->Get_FullConstrainParticles(), locStepVertexInfo->Get_OnlyConstrainTimeParticles(), locStepVertexInfo->Get_NoConstrainParticles(), locX4));
			else
				locSortedVertexConstraints.push_back(Make_VertexConstraint(locStepVertexInfo->Get_FullConstrainParticles(), locStepVertexInfo->Get_NoConstrainParticles(), locX4.Vect()));
		}
	}
	locAllConstraints.insert(locSortedVertexConstraints.begin(), locSortedVertexConstraints.end());

	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: All Constraints Created." << endl;

	return locAllConstraints;
}

/************************************************************** CONSTRAINT PREDICTORS **************************************************************/

//These functions are necessary to determine:
	//Whether each neutral is in a vertex constraint or not (may not be enough particles to constrain that particular vertex)
	//The pull terms needed when creating histograms
	//The constraint strings for the confidence level histogram

set<pair<int, int> > DKinFitUtils_GlueX::Get_KinFitVertexParticles(const DReaction* locReaction) const
{
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	if((locKinFitType == d_NoFit) || (locKinFitType == d_P4Fit))
		return set<pair<int, int> >();

	//predict vertices
	string locDummyString;
	size_t locNumConstraints = 0;
pair<size_t, string> DKinFitUtils_GlueX::Predict_VertexConstraints(const DReactionVertexInfo* locReactionVertexInfo, bool locSpacetimeFitFlag)
	locVertices = Predict_VertexConstraints(locReaction, locVertices, false, locNumConstraints, locDummyString); //false: doesn't matter: not used

	//merge vertices together
	set<pair<int, int> > locVertexParticles;
	for(size_t loc_i = 0; loc_i < locVertices.size(); ++loc_i)
		locVertexParticles.insert(locVertices[loc_i].begin(), locVertices[loc_i].end());

	return locVertexParticles;
}

string DKinFitUtils_GlueX::Get_ConstraintInfo(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, size_t& locNumConstraints, size_t& locNumUnknowns) const
{
	//returns constraint string, sets # constraints & unknowns
	//ASSUMES: Detected particles have non-zero errors!
	string locAllConstraintsString;
	locNumConstraints = 0;
	locNumUnknowns = 0;

	//P4 & Mass Constraints
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//Mass Constraints
		//loop over steps, but skip init step (if open-ended decaying, do p4 constraint instead)
		map<size_t, string> locMassConstraintStrings;
		for(size_t loc_i = 1; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
		{
			const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
			if(!locReactionStep->Get_KinFitConstrainInitMassFlag())
				continue; //don't apply mass constraint to this step

			Particle_t locPID = locReactionStep->Get_InitialParticleID();
			if(!IsFixedMass(locPID))
				continue; //don't apply mass constraint to this step

			locMassConstraintStrings[loc_i] = string("#it{m}_{") + string(ParticleName_ROOT(locPID)) + string("}");
		}

		//P4 Constraint //don't do if inclusive reaction
		Particle_t locMissingPID = Unknown;
		bool locMissingParticleUsedFlag = locReaction->Get_MissingPID(locMissingPID);
		if(!locMissingParticleUsedFlag || (locMissingPID != Unknown))
		{
			//OK, now, check to see if the system is overly constrained: 
				//there must be at least one particle with non-zero errors in the p4 constraint that is NOT in a mass constraint
			bool locNonZeroErrorFlag = false;

			//Find step used for p4 constraint: the one with missing/open-ended-decaying particle (if any: else first step)
			int locDefinedParticleStepIndex = locReaction->Get_DefinedParticleStepIndex();
			int locP4StepIndex = (locDefinedParticleStepIndex >= 0) ? locDefinedParticleStepIndex : 0;
			const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locP4StepIndex);

			set<size_t> locP4ConstrainedParticleSteps;

			//check if initial & final states if have non-zero cov errors:
			if((locP4StepIndex == 0) && (locReactionStep->Get_TargetParticleID() != Unknown) && dWillBeamHaveErrorsFlag)
				locNonZeroErrorFlag = true; //beam: we're good
			else if((locP4StepIndex != 0) && (locMassConstraintStrings.find(locP4StepIndex) == locMassConstraintStrings.end()))
				locNonZeroErrorFlag = true; //decaying particle, but mass not constrained: we're good (unless it's e.g. an omega. ugh.)
			else //check final state
			{
				//check if final state has non-zero cov errors (detected):
				for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalParticleIDs(); ++loc_i)
				{
					if(locReactionStep->Get_MissingParticleIndex() == int(loc_i))
						continue; //missing

					int locDecayStepIndex = locReaction->Get_DecayStepIndex(locP4StepIndex, loc_i);
					if(locDecayStepIndex > 0) //decaying
					{
						if(locMassConstraintStrings.find(locDecayStepIndex) != locMassConstraintStrings.end())
						{
							locP4ConstrainedParticleSteps.insert(locDecayStepIndex);
							continue; //mass constrained
						}
					}

					locNonZeroErrorFlag = true;
					break; //either detected particle, or unconstrained decaying particle: we're good (unless it's e.g. an omega. ugh.)
				}
			}

			bool locIncludeP4ConstraintFlag = true;
			if(!locNonZeroErrorFlag)
			{
				//system is over-constrained: we must delete a constraint
				//if there is a missing/open-ended particle: delete a mass constraint; else, delete the p4 constraint
				if((locDefinedParticleStepIndex < 0) || locP4ConstrainedParticleSteps.empty())
					locIncludeP4ConstraintFlag = false; //remove the p4 constraint
				else //remove a mass constraint: delete the one in the earliest step (so consistent) (will be by missing mass if present)
				{
					size_t locEarliestStepIndex = *locP4ConstrainedParticleSteps.begin();
					locMassConstraintStrings.erase(locEarliestStepIndex); //remove constraint
				}
			}

			//Finally, add p4 constraint
			if(locIncludeP4ConstraintFlag)
			{
				locAllConstraintsString = "#it{p}^{4}";
				locNumConstraints += 4;
				if(locDefinedParticleStepIndex >= 0)
					locNumUnknowns += 3;
			}
		}

		//Finally, add remaining mass constraints
		locNumConstraints += locMassConstraintStrings.size();
		map<size_t, string>::iterator locStringIterator = locMassConstraintStrings.begin();
		for(; locStringIterator != locMassConstraintStrings.end(); ++locStringIterator)
		{
			if(locAllConstraintsString != "")
				locAllConstraintsString += ", ";
			locAllConstraintsString += locStringIterator->second;
		}
	}

	//Create Vertex Constraints
	//VERTEX OR SPACETIME: Group particles by detached vertex (one deque for each constraint/vertex)
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		bool locSpacetimeFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));
		auto locConstraintPair = Predict_VertexConstraints(locReactionVertexInfo, locSpacetimeFitFlag);
		if(locConstraintPair.first > 0)
		{
			if(locAllConstraintsString != "")
				locAllConstraintsString += ", ";
			locAllConstraintsString += locConstraintPair.second;

			locNumConstraints += 2*DAnalysis::Get_FullConstrainParticles(locReactionVertexInfo).size();
			locNumConstraints += DAnalysis::Get_OnlyConstrainTimeParticles(locReactionVertexInfo).size();
			if(locSpacetimeFitFlag)
				locNumUnknowns += 4*locConstraintPair.first;
			else
				locNumUnknowns += 3*locConstraintPair.first;
		}
	}

	return locAllConstraintsString;
}

pair<size_t, string> DKinFitUtils_GlueX::Predict_VertexConstraints(const DReactionVertexInfo* locReactionVertexInfo, bool locSpacetimeFitFlag) const
{
	//returned: #constraints, constraint string
	size_t locNumConstraints = 0;
	string locAllConstraintString;

	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	for(auto& locVertexInfo : locStepVertexInfos)
	{
		if(locVertexInfo->Get_DanglingVertexFlag())
			continue;

		auto locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles();
		if(locSpacetimeFitFlag)
		{
			locNumConstraints += 3*locFullConstrainParticles.size();
			locNumConstraints += locVertexInfo->Get_OnlyConstrainTimeParticles().size();
		}
		else //vertex only
			locNumConstraints += 2*locFullConstrainParticles.size();

		//adjust if beamline not included in vertex fit
		if(!Get_IncludeBeamlineInVertexFitFlag() && locVertexInfo->Get_ProductionVertexFlag())
			locNumConstraints -= 2*locVertexInfo->Get_FullConstrainParticles(d_InitialState).size();

		//add to the full constraint string
		if(locAllConstraintString != "")
			locAllConstraintString += ", ";
		locAllConstraintString += Build_VertexConstraintString(locVertexInfo, locSpacetimeFitFlag);
	}

	return make_pair(locNumConstraints, locAllConstraintString);
}

string DKinFitUtils_GlueX::Build_VertexConstraintString(const DReactionStepVertexInfo* locVertexInfo, bool locSpacetimeFitFlag)
{
	auto locReaction = locVertexInfo->Get_Reaction();
	string locConstraintString = locSpacetimeFitFlag ? "#it{x}^{4}_{" : "#it{x}^{3}_{";

	//initial state
	auto locParticles = locVertexInfo->Get_Particles(d_InitialState);
	auto locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles(d_InitialState);
	for(auto locIndices : locParticles)
	{
		auto locStep = locReaction->Get_ReactionStep(locIndices.first);
		Particle_t locPID = locStep->Get_PID(locIndices.second);
		string locParticleString = ParticleName_ROOT(locPID);
		if(locIndices.second == locStep->Get_MissingParticleIndex())
			locConstraintString += string("(") + locParticleString + string(")"); //missing
		else if(std::binary_search(locFullConstrainParticles.begin(), locFullConstrainParticles.end(), locIndices)) //constraining
		{
			//adjust if beamline not included in vertex fit
			if(!Get_IncludeBeamlineInVertexFitFlag() && locVertexInfo->Get_ProductionVertexFlag())
				locConstraintString += locParticleString; //no-constrain
			else //constrain
				locConstraintString += string("#color[4]{") + locParticleString + string("}"); //blue
		}
		else //no-constrain
			locConstraintString += locParticleString; //plain
	}

	//final state
	locConstraintString += "#rightarrow";
	locParticles = locVertexInfo->Get_Particles(d_FinalState);
	locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles(d_FinalState);
	auto locOnlyConstrainTimeParticles = locVertexInfo->Get_OnlyConstrainTimeParticles();
	for(auto locIndices : locParticles)
	{
		auto locStep = locReaction->Get_ReactionStep(locIndices.first);
		Particle_t locPID = locStep->Get_PID(locIndices.second);
		string locParticleString = ParticleName_ROOT(locPID);

		if(locIndices.second == locStep->Get_MissingParticleIndex())
			locConstraintString += string("(") + locParticleString + string(")"); //missing
		else if(std::binary_search(locFullConstrainParticles.begin(), locFullConstrainParticles.end(), locIndices)) //constraining
			locConstraintString += string("#color[4]{") + locParticleString + string("}"); //blue
		else if(locSpacetimeFitFlag && std::binary_search(locOnlyConstrainTimeParticles.begin(), locOnlyConstrainTimeParticles.end(), locIndices)) //time-only
			locConstraintString += string("#color[3]{") + locParticleString + string("}"); //green
		else //no-constrain
			locConstraintString += locParticleString; //plain
	}

	locConstraintString += string("}"); //end of constraint

	return locConstraintString;
}

/*************************************************************** CALCULATION ROUTINES **************************************************************/

bool DKinFitUtils_GlueX::Propagate_TrackInfoToCommonVertex(DKinematicData* locKinematicData, DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi)
{
	//locKinematicData must be generated from locKinFitParticle
		//this function should only be used on decaying particles involved in two vertex fits:
			//propagates the track information from the vertex at which it is DEFINED to the OTHER vertex (e.g. production -> decay)

	TVector3 locMomentum;
	TLorentzVector locSpacetimeVertex;
	pair<double, double> locPathLengthPair;
	TMatrixFSym* locCovarianceMatrix = Get_SymMatrixResource(7);
	if(!DKinFitUtils::Propagate_TrackInfoToCommonVertex(locKinFitParticle, locVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, locCovarianceMatrix))
		return false;

	locKinematicData->setMomentum(DVector3(locMomentum.X(),locMomentum.Y(),locMomentum.Z()));
	locKinematicData->setPosition(DVector3(locSpacetimeVertex.Vect().X(),locSpacetimeVertex.Vect().Y(),locSpacetimeVertex.Vect().Z()));
	locKinematicData->setTime(locSpacetimeVertex.T());
	locKinematicData->setErrorMatrix(locCovarianceMatrix);
	locKinematicData->setPathLength(locPathLengthPair.first, locPathLengthPair.second);
	return true;
}
