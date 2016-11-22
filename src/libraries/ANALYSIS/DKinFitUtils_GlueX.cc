#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DKinFitUtils_GlueX.h"

/******************************************************************** INITIALIZE *******************************************************************/

DKinFitUtils_GlueX::DKinFitUtils_GlueX(const DMagneticFieldMap* locMagneticFieldMap, const DAnalysisUtilities* locAnalysisUtilities) : 
dMagneticFieldMap(locMagneticFieldMap), dAnalysisUtilities(locAnalysisUtilities)
{
	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);
	dWillBeamHaveErrorsFlag = false; //Until fixed!
}

DKinFitUtils_GlueX::DKinFitUtils_GlueX(JEventLoop* locEventLoop)
{
	locEventLoop->GetSingle(dAnalysisUtilities);

	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	dMagneticFieldMap = locApplication->GetBfield(locEventLoop->GetJEvent().GetRunNumber());

	gPARMS->SetDefaultParameter("KINFIT:LINKVERTICES", dLinkVerticesFlag);
	dWillBeamHaveErrorsFlag = false; //Until fixed!
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

	Set_MaxMatrixDSymPoolSize(10*locNumReactions*locExpectedNumCombos*2*5); //extra x5: to be safe
}

/*********************************************************** OVERRIDE BASE CLASS FUNCTIONS *********************************************************/

void DKinFitUtils_GlueX::Reset_NewEvent(void)
{
	dParticleMap_SourceToInput_Beam.clear();
	dParticleMap_SourceToInput_DetectedParticle.clear();
	dParticleMap_SourceToInput_DetectedParticleFromDecay.clear();
	dParticleMap_SourceToInput_Shower.clear();
	dParticleMap_SourceToInput_Target.clear();
	dParticleMap_SourceToInput_Decaying.clear();
	dParticleMap_SourceToInput_Missing.clear();

	dParticleMap_InputToSource_JObject.clear();
	dParticleMap_InputToSource_Decaying.clear();

	DKinFitUtils::Reset_NewEvent();
}

bool DKinFitUtils_GlueX::Get_IncludeBeamlineInVertexFitFlag(void) const
{
	return false; //at least until covariance matrix is set for beam photons
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
		locSpacetimeVertex, locMomentum, &(locBeamPhoton->errorMatrix()));
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
	TMatrixDSym locCovarianceMatrix = locBeamPhoton->errorMatrix();
	locCovarianceMatrix(6, 6) = locEventRFBunch->dTimeVariance;
	//zero the correlation terms
	for(int loc_i = 0; loc_i < 6; ++loc_i)
	{
		locCovarianceMatrix(6, loc_i) = 0.0;
		locCovarianceMatrix(loc_i, 6) = 0.0;
	}

	DKinFitParticle* locKinFitParticle = DKinFitUtils::Make_BeamParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), 
		locSpacetimeVertex, locMomentum, &locCovarianceMatrix);
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
		locSpacetimeVertex, locMomentum, &(locKinematicData->errorMatrix()));
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

	if(dParticleMap_SourceToInput_DetectedParticleFromDecay.find(locDecayingKinFitParticle) != dParticleMap_SourceToInput_DetectedParticleFromDecay.end())
		return dParticleMap_SourceToInput_DetectedParticleFromDecay[locDecayingKinFitParticle]; //not unique, return existing

	DKinFitParticle* locDetectedKinFitParticle = DKinFitUtils::Make_DetectedParticle(locDecayingKinFitParticle->Get_PID(), 
		locDecayingKinFitParticle->Get_Charge(), locDecayingKinFitParticle->Get_Mass(), locDecayingKinFitParticle->Get_SpacetimeVertex(), 
		locDecayingKinFitParticle->Get_Momentum(), locDecayingKinFitParticle->Get_CovarianceMatrix());

	dParticleMap_SourceToInput_DetectedParticleFromDecay[locDecayingKinFitParticle] = locDetectedKinFitParticle;
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
		locNeutralShower->dEnergy, &(locNeutralShower->dCovarianceMatrix));

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

set<DKinFitConstraint*> DKinFitUtils_GlueX::Create_Constraints(const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain, DKinFitType locKinFitType, deque<DKinFitConstraint_Vertex*>& locSortedVertexConstraints)
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
	if((locKinFitType == d_VertexFit) || (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		bool locSpacetimeFitFlag = ((locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit));
		locSortedVertexConstraints = Create_VertexConstraints(locKinFitChain, locSpacetimeFitFlag);
	}
	locAllConstraints.insert(locSortedVertexConstraints.begin(), locSortedVertexConstraints.end());

	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: All Constraints Created." << endl;

	return locAllConstraints;
}

/********************************************************* MAKE INITIAL SPACETIME GUESSES **********************************************************/

void DKinFitUtils_GlueX::Set_SpacetimeGuesses(const deque<DKinFitConstraint_Vertex*>& locSortedVertexConstraints, bool locIsP4FitFlag)
{
	//loop through vertices, determining initial guesses
	map<DKinFitParticle*, DKinFitParticle*> locDecayingToDetectedParticleMap; //input decaying particle -> new detected particle
	for(size_t loc_i = 0; loc_i < locSortedVertexConstraints.size(); ++loc_i)
	{
		//get constraint
		DKinFitConstraint_Vertex* locOrigVertexConstraint = locSortedVertexConstraints[loc_i];
		DKinFitConstraint_Spacetime* locOrigSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOrigVertexConstraint);

		DKinFitConstraint_Vertex* locActiveVertexConstraint = locOrigVertexConstraint;
		DKinFitConstraint_Spacetime* locActiveSpacetimeConstraint = locOrigSpacetimeConstraint;

		/**************************************************** SUBSTITUTE FOR DECAYING PARTICLES ******************************************************/

		//If a decaying particle was previously reconstructed, substitute for it so can do a stand-alone vertex fit
		bool locAttemptFitFlag = true; //if set to false in Build_NewConstraint(), will not attempt fit (a previous fit failed)
		set<DKinFitParticle*> locFullConstrainSet = locOrigVertexConstraint->Get_FullConstrainParticles();
		size_t locNumDecayingConstrainParticles = 0;
		set<DKinFitParticle*>::iterator locParticleIterator = locFullConstrainSet.begin();
		for(; locParticleIterator != locFullConstrainSet.end(); ++locParticleIterator)
		{
			if((*locParticleIterator)->Get_KinFitParticleType() == d_DecayingParticle)
				++locNumDecayingConstrainParticles;
		}
		if(locNumDecayingConstrainParticles > 0)
		{
			//true: skip bad decaying particles (those whose reconstruction-fits failed) if at all possible
				//if cannot skip, will set locAttemptFitFlag to false
			locActiveVertexConstraint = Build_NewConstraint(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locAttemptFitFlag, true);
			locActiveSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locActiveVertexConstraint);
		}

		/*********************************************************** INITIAL VERTEX GUESS ************************************************************/

		//do crude vertex guess: point on DOCA-line between the two closest (by doca) particles

		//get particles
		locFullConstrainSet = locActiveVertexConstraint->Get_FullConstrainParticles();
		deque<DKinFitParticle*> locFullConstrainDeque;
		std::copy(locFullConstrainSet.begin(), locFullConstrainSet.end(), std::back_inserter(locFullConstrainDeque));

		//get guess
		DVector3 locDVertexGuess = dAnalysisUtilities->Calc_CrudeVertex(locFullConstrainDeque);
		TVector3 locVertexGuess(locDVertexGuess.X(), locDVertexGuess.Y(), locDVertexGuess.Z());
		if(dDebugLevel > 20)
			cout << "init vertex guess = " << locVertexGuess.X() << ", " << locVertexGuess.Y() << ", " << locVertexGuess.Z() << endl;

		//set guess
		locOrigVertexConstraint->Set_InitVertexGuess(locVertexGuess);
		locActiveVertexConstraint->Set_InitVertexGuess(locVertexGuess);

		/************************************************************ INITIAL TIME GUESS *************************************************************/

		double locTimeGuess = 0.0;
		if(locActiveSpacetimeConstraint != NULL)
		{
			locTimeGuess = Calc_TimeGuess(locActiveSpacetimeConstraint, locDVertexGuess);
			locOrigSpacetimeConstraint->Set_InitTimeGuess(locTimeGuess);
			locActiveSpacetimeConstraint->Set_InitTimeGuess(locTimeGuess);
		}

		/************************************************************** VERTEX/TIME FIT **************************************************************/

		if((locSortedVertexConstraints.size() == 1) && !locIsP4FitFlag)
			return; //only one constraint: don't do primary fit here

		//Check if should attempt fit
		if(!locAttemptFitFlag)
		{
			//No: The results from a needed previous fit failed. Create new decaying particles and move on.
			TLorentzVector locSpacetimeVertex(locVertexGuess, locTimeGuess);
			Construct_DetectedDecayingParticle_NoFit(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locSpacetimeVertex);
			continue;
		}

		//Fit & Update vertex guess
		dKinFitter->Reset_NewFit();
		dKinFitter->Add_Constraint(locActiveVertexConstraint);
		bool locFitStatus = dKinFitter->Fit_Reaction();

		/************************************************************ UPDATE WITH RESULTS ************************************************************/

		if(locFitStatus) //True if fit succeeded
		{
			//Update vertex guess
			DKinFitConstraint_Vertex* locNewVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*dKinFitter->Get_KinFitConstraints().begin());
			locOrigVertexConstraint->Set_InitVertexGuess(locNewVertexConstraint->Get_CommonVertex());

			//Update time guess (if time fit)
			if(locOrigSpacetimeConstraint != NULL)
			{
				DKinFitConstraint_Spacetime* locNewSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locNewVertexConstraint);
				locOrigSpacetimeConstraint->Set_InitTimeGuess(locNewSpacetimeConstraint->Get_CommonTime());
			}

			//create detected particles out of reconstructed decaying particles in this constraint
			set<DKinFitParticle*> locOutputKinFitParticles = dKinFitter->Get_KinFitParticles();
			set<DKinFitParticle*>::iterator locResultIterator = locOutputKinFitParticles.begin();
			for(; locResultIterator != locOutputKinFitParticles.end(); ++locResultIterator)
			{
				if((*locResultIterator)->Get_KinFitParticleType() != d_DecayingParticle)
					continue;

				set<DKinFitParticle*> locAllVertexParticles = locNewVertexConstraint->Get_AllParticles();
				if(locAllVertexParticles.find(*locResultIterator) == locAllVertexParticles.end())
					continue; //not used in this constraint: vertex not yet defined

				DKinFitParticle* locInputKinFitParticle = Get_InputKinFitParticle(*locResultIterator);
				locDecayingToDetectedParticleMap[locInputKinFitParticle] = Make_DetectedParticle(*locResultIterator);
			}
		}
		else //fit failed, but still need to reconstruct decaying particles for next vertex-find step
		{
			TLorentzVector locSpacetimeVertex(locVertexGuess, locTimeGuess);
			Construct_DetectedDecayingParticle_NoFit(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locSpacetimeVertex);
		}
	}
}

void DKinFitUtils_GlueX::Construct_DetectedDecayingParticle_NoFit(DKinFitConstraint_Vertex* locOrigVertexConstraint, map<DKinFitParticle*, DKinFitParticle*>& locDecayingToDetectedParticleMap, TLorentzVector locSpacetimeVertexGuess)
{
	//get "decaying" no-constrain decaying particles
	set<DKinFitParticle*> locNoConstrainParticles = locOrigVertexConstraint->Get_NoConstrainParticles();
	set<DKinFitParticle*>::iterator locNoConstrainIterator = locNoConstrainParticles.begin();
	for(; locNoConstrainIterator != locNoConstrainParticles.end(); ++locNoConstrainIterator)
	{
		DKinFitParticle* locInputKinFitParticle = *locNoConstrainIterator;
		if(locInputKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //not a decaying particle

		//create a new one
		TLorentzVector locP4 = Calc_DecayingP4_ByP3Derived(locInputKinFitParticle, true, true);
		TMatrixDSym locCovarianceMatrix(7);
		locCovarianceMatrix(0, 0) = -1.0; //signal that you shouldn't do fits that need this particle
		DKinFitParticle* locDetectedKinFitParticle = Make_DetectedParticle(locInputKinFitParticle->Get_PID(), 
			locInputKinFitParticle->Get_Charge(), locInputKinFitParticle->Get_Mass(), locSpacetimeVertexGuess, locP4.Vect(), &locCovarianceMatrix);

		//register it
		locDecayingToDetectedParticleMap[locInputKinFitParticle] = locDetectedKinFitParticle;
	}
}

DKinFitConstraint_Vertex* DKinFitUtils_GlueX::Build_NewConstraint(DKinFitConstraint_Vertex* locOrigVertexConstraint, const map<DKinFitParticle*, DKinFitParticle*>& locDecayingToDetectedParticleMap, bool& locAttemptFitFlag, bool locSkipBadDecayingFlag)
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX::Build_NewConstraint()" << endl;
	set<DKinFitParticle*> locNewDetectedParticles, locUsedDecayingParticles;

	//get "detected" versions of reconstructed decaying particles
	set<DKinFitParticle*> locFullConstrainParticles = locOrigVertexConstraint->Get_FullConstrainParticles();
	set<DKinFitParticle*>::iterator locFullConstrainIterator = locFullConstrainParticles.begin();
	set<DKinFitParticle*> locNewFullConstrainParticles;
	for(; locFullConstrainIterator != locFullConstrainParticles.end(); ++locFullConstrainIterator)
	{
		DKinFitParticle* locInputKinFitParticle = *locFullConstrainIterator;
		if(locInputKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		{
			locNewFullConstrainParticles.insert(locInputKinFitParticle);
			continue; //not a decaying particle
		}

		map<DKinFitParticle*, DKinFitParticle*>::const_iterator locDecayIterator = locDecayingToDetectedParticleMap.find(locInputKinFitParticle);
		if(locDecayIterator == locDecayingToDetectedParticleMap.end())
		{
			if(!locSkipBadDecayingFlag)
				locAttemptFitFlag = false; //cannot fit. however, still get initial guess
			continue; //try to see if can fit without this particle
		}
		DKinFitParticle* locDetectedDecayingParticle = locDecayIterator->second;
		if((*(locDetectedDecayingParticle->Get_CovarianceMatrix()))(0, 0) < 0.0)
		{
			if(locSkipBadDecayingFlag)
				continue; //try to see if can fit without this particle
			else
				locAttemptFitFlag = false; //cannot fit. however, still get initial guess
		}

		//reconstructed particle found
		locNewFullConstrainParticles.insert(locDetectedDecayingParticle);
	}

	//Check if have enough particles
	if(locNewFullConstrainParticles.size() < 2) //cannot fit: try again, using non-fit decaying particles
		return Build_NewConstraint(locOrigVertexConstraint, locDecayingToDetectedParticleMap, locAttemptFitFlag, false);

	//create new constraint, this time with the new detected particles
	set<DKinFitParticle*> locNoConstrainParticles = locOrigVertexConstraint->Get_NoConstrainParticles();
	DKinFitConstraint_Spacetime* locOrigSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOrigVertexConstraint);
	if(locOrigSpacetimeConstraint == NULL) //vertex fit
		return Make_VertexConstraint(locNewFullConstrainParticles, locNoConstrainParticles);
	else
		return Make_SpacetimeConstraint(locNewFullConstrainParticles, locOrigSpacetimeConstraint->Get_OnlyConstrainTimeParticles(), locNoConstrainParticles);
}

double DKinFitUtils_GlueX::Calc_TimeGuess(const DKinFitConstraint_Spacetime* locConstraint, DVector3 locVertexGuess)
{
	set<DKinFitParticle*> locTimeFindParticles = locConstraint->Get_FullConstrainParticles();
	set<DKinFitParticle*> locOnlyTimeFindParticles = locConstraint->Get_OnlyConstrainTimeParticles();

	//see if can find the beam particle: if so, use the rf time
	set<DKinFitParticle*>& locSearchBeamParticles = Get_IncludeBeamlineInVertexFitFlag() ? locTimeFindParticles : locOnlyTimeFindParticles;
	set<DKinFitParticle*>::iterator locParticleIterator = locSearchBeamParticles.begin();
	for(; locParticleIterator != locSearchBeamParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_BeamParticle)
			continue;

		//have the beam particle: use the rf time (propagate it to the vertex)
		double locDeltaZ = locVertexGuess.Z() - locKinFitParticle->Get_Position().Z();
		return (locKinFitParticle->Get_Time() + locDeltaZ/29.9792458);
	}

	//propagate each track time to the DOCA to the init vertex guess and average them
	locTimeFindParticles.insert(locOnlyTimeFindParticles.begin(), locOnlyTimeFindParticles.end());

	//build vector
	deque<DKinFitParticle*> locTimeFindParticleVector;
	std::copy(locTimeFindParticles.begin(), locTimeFindParticles.end(), std::back_inserter(locTimeFindParticleVector));

	return dAnalysisUtilities->Calc_CrudeTime(locTimeFindParticleVector, locVertexGuess);
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
	deque<set<pair<int, int> > > locVertices = Setup_VertexPredictions(locReaction);
	string locDummyString;
	size_t locNumConstraints = 0;
	locVertices = Predict_VertexConstraints(locReaction, locVertices, false, locNumConstraints, locDummyString); //false: doesn't matter: not used

	//merge vertices together
	set<pair<int, int> > locVertexParticles;
	for(size_t loc_i = 0; loc_i < locVertices.size(); ++loc_i)
		locVertexParticles.insert(locVertices[loc_i].begin(), locVertices[loc_i].end());

	return locVertexParticles;
}

string DKinFitUtils_GlueX::Get_ConstraintInfo(const DReaction* locReaction, DKinFitType locKinFitType, size_t& locNumConstraints, size_t& locNumUnknowns) const
{
	//returns constraint string, sets # constraints & unknowns
	//ASSUMES: Detected particles have non-zero errors!
	string locAllConstraintsString;
	locNumConstraints = 0;
	locNumUnknowns = 0;

	//P4 & Mass Constraints
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
		deque<set<pair<int, int> > > locVertices = Setup_VertexPredictions(locReaction);

		string locVertexConstraintString;
		size_t locNumVertexConstraints = 0;
		locVertices = Predict_VertexConstraints(locReaction, locVertices, locSpacetimeFitFlag, locNumVertexConstraints, locVertexConstraintString);

		if(!locVertices.empty())
		{
			if(locAllConstraintsString != "")
				locAllConstraintsString += ", ";
			locAllConstraintsString += locVertexConstraintString;

			locNumConstraints += locNumVertexConstraints;
			if(locSpacetimeFitFlag)
				locNumUnknowns += 4*locVertices.size();
			else
				locNumUnknowns += 3*locVertices.size();
		}
	}

	return locAllConstraintsString;
}

deque<set<pair<int, int> > > DKinFitUtils_GlueX::Setup_VertexPredictions(const DReaction* locReaction) const
{
	//create decay maps
	map<pair<int, int>, int> locDecayMap_ParticleToDecayStep;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(loc_i);
		for(size_t loc_j = 0; loc_j < locReactionStep->Get_NumFinalParticleIDs(); ++loc_j)
		{
			int locDecayStepIndex = locReaction->Get_DecayStepIndex(loc_i, loc_j);
			if(locDecayStepIndex < 0)
				continue;
			locDecayMap_ParticleToDecayStep[pair<int, int>(loc_i, loc_j)] = locDecayStepIndex; //store step where this particle decays
		}
	}

	//no choice but to repeat what's been done, but without DKinFitChain:
	deque<set<pair<int, int> > > locAllVertices; //one pair for each vertex: particles, constraint string
	set<size_t> locIncludedStepIndices;

	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		if(locIncludedStepIndices.find(loc_i) != locIncludedStepIndices.end())
			continue; //already did this step

		//Start a new vertex, and save when done
		set<pair<int, int> > locVertexParticles;
		Setup_VertexPrediction(locReaction, loc_i, locVertexParticles, locDecayMap_ParticleToDecayStep, locIncludedStepIndices);
		locAllVertices.push_back(locVertexParticles);
	}

	return locAllVertices;
}

void DKinFitUtils_GlueX::Setup_VertexPrediction(const DReaction* locReaction, size_t locStepIndex, set<pair<int, int> >& locVertexParticles, const map<pair<int, int>, int>& locDecayMap_ParticleToDecayStep, set<size_t>& locIncludedStepIndices) const
{
	bool locStartNewVertexFlag = locVertexParticles.empty();
	locIncludedStepIndices.insert(locStepIndex);

	//get the step
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);

	//if new constraint: loop over the initial particles: add to constraint string
	if(locStartNewVertexFlag)
	{
		locVertexParticles.insert(pair<int, int>(locStepIndex, -2)); //beam/decaying particle
		if(locReactionStep->Get_TargetParticleID() != Unknown)
			locVertexParticles.insert(pair<int, int>(locStepIndex, -1)); //target
	}

	//loop over final particles: add to the vertex constraint, dive through decaying particles that decay in-place
		//if decaying in-place: don't add (would add to constraint, but not here (purely internal))
	for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalParticleIDs(); ++loc_i)
	{
		pair<int, int> locParticlePair(locStepIndex, loc_i);
		Particle_t locPID = locReactionStep->Get_FinalParticleID(loc_i);

		//check if particle decays, and if so, if in-place
		map<pair<int, int>, int>::const_iterator locDecayIterator = locDecayMap_ParticleToDecayStep.find(locParticlePair);
		int locDecayStepIndex = (locDecayIterator == locDecayMap_ParticleToDecayStep.end()) ? -1 : locDecayIterator->second;
		if((locDecayStepIndex >= 0) && !IsDetachedVertex(locPID))
		{
			//yes: combine with the decay products
			Setup_VertexPrediction(locReaction, locDecayStepIndex, locVertexParticles, locDecayMap_ParticleToDecayStep, locIncludedStepIndices);
			continue;
		}
		else //does not decay, or at least, not in-place
			locVertexParticles.insert(locParticlePair);
	}
}

deque<set<pair<int, int> > > DKinFitUtils_GlueX::Predict_VertexConstraints(const DReaction* locReaction, deque<set<pair<int, int> > > locAllVertices, bool locSpacetimeFitFlag, size_t& locNumConstraints, string& locAllConstraintString) const
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: Create vertex constraints." << endl;

	locNumConstraints = 0;

	//resolve links between vertex & time fits (decaying particles), create constraints, sort them, and return them
		//sort: the order in which they are defined (as required by the decaying particles they're using)
	deque<set<pair<int, int> > > locAllVertexParticles;

	//initialize particle groupings
	deque<set<pair<int, int> > > locAllFullConstrainParticles, locAllDecayingParticles, locAllOnlyConstrainTimeParticles, locAllNoConstrainParticles;
	for(size_t loc_i = 0; loc_i < locAllVertices.size(); ++loc_i)
	{
		set<pair<int, int> > locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles;
		Group_VertexParticles(locReaction, locSpacetimeFitFlag, locAllVertices[loc_i], locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles);
		locAllFullConstrainParticles.push_back(locFullConstrainParticles);
		locAllDecayingParticles.push_back(locDecayingParticles);
		locAllOnlyConstrainTimeParticles.push_back(locOnlyConstrainTimeParticles);
		locAllNoConstrainParticles.push_back(locNoConstrainParticles);
	}

	//loop over vertex-constraints-to-sort:
		//find which constraints decaying particles should be defined-by/constrained-to
		//find order in which constraints need to be constrained
	size_t locConstraintIndex = 0;
	bool locProgessMadeFlag = false;
	set<pair<int, int> > locDefinedDecayingParticles;
	while(!locAllFullConstrainParticles.empty())
	{
		if(locConstraintIndex == locAllFullConstrainParticles.size())
		{
			//made a full loop through
			if(!locProgessMadeFlag)
				break; //no progress made: cannot constrain remaining vertices
			//reset for next pass through
			locConstraintIndex = 0;
			locProgessMadeFlag = false;
			continue;
		}

		//find which decaying particles at this vertex have been previously defined
		set<pair<int, int> > locVertexDecayingParticles_Defined;
		if(dLinkVerticesFlag)
			set_intersection(locAllDecayingParticles[locConstraintIndex].begin(), locAllDecayingParticles[locConstraintIndex].end(),
				locDefinedDecayingParticles.begin(), locDefinedDecayingParticles.end(),
				inserter(locVertexDecayingParticles_Defined, locVertexDecayingParticles_Defined.begin()));

		//see if enough defined particles to constrain vertex
		if(locVertexDecayingParticles_Defined.size() + locAllFullConstrainParticles[locConstraintIndex].size() < 2)
		{
			++locConstraintIndex;
			continue; //nope
		}

		//find which decaying particles at this vertex have NOT been previously defined
		set<pair<int, int> > locVertexDecayingParticles_NotDefined;
		set_difference(locAllDecayingParticles[locConstraintIndex].begin(), locAllDecayingParticles[locConstraintIndex].end(),
			locDefinedDecayingParticles.begin(), locDefinedDecayingParticles.end(),
			inserter(locVertexDecayingParticles_NotDefined, locVertexDecayingParticles_NotDefined.begin()));

		//Add decaying particles to appropriate sets
		locAllFullConstrainParticles[locConstraintIndex].insert(locVertexDecayingParticles_Defined.begin(), locVertexDecayingParticles_Defined.end());
		locAllNoConstrainParticles[locConstraintIndex].insert(locVertexDecayingParticles_NotDefined.begin(), locVertexDecayingParticles_NotDefined.end());

		//Create vertex/spacetime constraint, with decaying particles as no-constrain particles
		set<pair<int, int> > locVertexParticles;
		if(locSpacetimeFitFlag)
		{
			locVertexParticles.insert(locAllOnlyConstrainTimeParticles[locConstraintIndex].begin(), locAllOnlyConstrainTimeParticles[locConstraintIndex].end());
			locNumConstraints += 3*locAllFullConstrainParticles[locConstraintIndex].size();
			locNumConstraints += locAllOnlyConstrainTimeParticles[locConstraintIndex].size();
		}
		else //vertex only
			locNumConstraints += 2*locAllFullConstrainParticles[locConstraintIndex].size();

		locVertexParticles.insert(locAllFullConstrainParticles[locConstraintIndex].begin(), locAllFullConstrainParticles[locConstraintIndex].end());
		locVertexParticles.insert(locAllNoConstrainParticles[locConstraintIndex].begin(), locAllNoConstrainParticles[locConstraintIndex].end());
		locAllVertexParticles.push_back(locVertexParticles);

		//The positions of these decaying particles are now defined: Can use to constrain vertices in later constraints
		if(dLinkVerticesFlag)
		{
			//since we need to match with particles in other constraints, save the OTHER index for the particle
				//if was in initial state, save final-state pair. and vice versa
			set<pair<int, int> >::iterator locDecayIterator = locVertexDecayingParticles_NotDefined.begin();
			for(; locDecayIterator != locVertexDecayingParticles_NotDefined.end(); ++locDecayIterator)
			{
				pair<int, int> locParticlePair = *locDecayIterator;
				if(locParticlePair.second < 0) //was in initial state: save final state
				{
					locParticlePair = locReaction->Get_InitialParticleDecayFromIndices(locParticlePair.first);
					locDefinedDecayingParticles.insert(locParticlePair);
				}
				else //was in final state: save final state
				{
					int locDecayStepIndex = locReaction->Get_DecayStepIndex(locParticlePair.first, locParticlePair.second);
					locDefinedDecayingParticles.insert(pair<int, int>(locDecayStepIndex, -2));
				}
			}
		}

		//add to the full constraint string
		if(locAllConstraintString != "")
			locAllConstraintString += ", ";
		locAllConstraintString += Build_VertexConstraintString(locReaction, locAllVertices[locConstraintIndex], locAllFullConstrainParticles[locConstraintIndex], locAllOnlyConstrainTimeParticles[locConstraintIndex], locAllNoConstrainParticles[locConstraintIndex], locSpacetimeFitFlag);

		//Erase this vertex from future consideration
		locAllVertices.erase(locAllVertices.begin() + locConstraintIndex);
		locAllFullConstrainParticles.erase(locAllFullConstrainParticles.begin() + locConstraintIndex);
		locAllDecayingParticles.erase(locAllDecayingParticles.begin() + locConstraintIndex);
		locAllOnlyConstrainTimeParticles.erase(locAllOnlyConstrainTimeParticles.begin() + locConstraintIndex);
		locAllNoConstrainParticles.erase(locAllNoConstrainParticles.begin() + locConstraintIndex);

		locProgessMadeFlag = true;
	}

	return locAllVertexParticles;
}

void DKinFitUtils_GlueX::Group_VertexParticles(const DReaction* locReaction, bool locSpacetimeFitFlag, const set<pair<int, int> >& locVertexParticles, set<pair<int, int> >& locFullConstrainParticles, set<pair<int, int> >& locDecayingParticles, set<pair<int, int> >& locOnlyConstrainTimeParticles, set<pair<int, int> >& locNoConstrainParticles) const
{
	//Decaying: Only those that can conceivably be used to constrain: All unless dLinkVerticesFlag disabled (then none)
	locFullConstrainParticles.clear();
	locDecayingParticles.clear();
	locOnlyConstrainTimeParticles.clear();
	locNoConstrainParticles.clear();

	set<pair<int, int> >::const_iterator locIterator = locVertexParticles.begin();
	for(; locIterator != locVertexParticles.end(); ++locIterator)
	{
		int locStepIndex = (*locIterator).first;
		int locParticleIndex = (*locIterator).second;

		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
		Particle_t locInitPID = locReactionStep->Get_InitialParticleID();
		int locDecayStepIndex = locReaction->Get_DecayStepIndex(locStepIndex, locParticleIndex);

		if(locParticleIndex == -1) //target
			locNoConstrainParticles.insert(*locIterator);
		else if((locParticleIndex == -2) && (locInitPID == Gamma)) //beam //ASSUMES BEAM TYPE
		{
			if(Get_IncludeBeamlineInVertexFitFlag())
				locFullConstrainParticles.insert(*locIterator);
			else
				locNoConstrainParticles.insert(*locIterator);
		}
		else if((locParticleIndex == -2) || (locDecayStepIndex >= 0)) //decaying
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.insert(*locIterator);
			else
				locNoConstrainParticles.insert(*locIterator);
		}
		else if(locParticleIndex == locReactionStep->Get_MissingParticleIndex()) //missing
			locNoConstrainParticles.insert(*locIterator);
		else if(ParticleCharge(locReactionStep->Get_FinalParticleID(locParticleIndex)) == 0) //detected neutral
		{
			if(locSpacetimeFitFlag)
				locOnlyConstrainTimeParticles.insert(*locIterator);
			else
				locNoConstrainParticles.insert(*locIterator);
		}
		else //detected charged
			locFullConstrainParticles.insert(*locIterator);
	}
}

string DKinFitUtils_GlueX::Build_VertexConstraintString(const DReaction* locReaction, const set<pair<int, int> >& locAllVertexParticles, set<pair<int, int> >& locFullConstrainParticles, set<pair<int, int> >& locOnlyConstrainTimeParticles, set<pair<int, int> >& locNoConstrainParticles, bool locSpacetimeFitFlag) const
{
	string locConstraintString;
	if(locSpacetimeFitFlag)
		locConstraintString += "#it{x}^{4}_{";
	else
		locConstraintString += "#it{x}^{3}_{";

	//initial particles
	bool locInitialStateFlag = true;
	set<pair<int, int> >::iterator locAllIterator = locAllVertexParticles.begin();
	for(; locAllIterator != locAllVertexParticles.end(); ++locAllIterator)
	{
		pair<int, int> locParticlePair = *locAllIterator;
		if(locInitialStateFlag && (locParticlePair.second >= 0))
		{
			//now on final state particles
			locInitialStateFlag = false;
			locConstraintString += "#rightarrow";
		}

		const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locParticlePair.first);
		Particle_t locPID = Unknown;
		if(locParticlePair.second == -2) //initial particle
			locPID = locReactionStep->Get_InitialParticleID();
		else if(locParticlePair.second == -1) //target particle
			locPID = locReactionStep->Get_TargetParticleID();
		else //final state
			locPID = locReactionStep->Get_FinalParticleID(locParticlePair.second);

		string locParticleString = ParticleName_ROOT(locPID);
		if((locParticlePair.second == locReactionStep->Get_MissingParticleIndex()) && (locParticlePair.second != -1)) //-1 = target
			locConstraintString += string("(") + locParticleString + string(")"); //missing
		else if(locFullConstrainParticles.find(locParticlePair) != locFullConstrainParticles.end()) //constraining
			locConstraintString += string("#color[4]{") + locParticleString + string("}"); //blue
		else if(locOnlyConstrainTimeParticles.find(locParticlePair) != locOnlyConstrainTimeParticles.end()) //time-only
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
	TMatrixDSym* locCovarianceMatrix = Get_MatrixDSymResource(7);
	if(!DKinFitUtils::Propagate_TrackInfoToCommonVertex(locKinFitParticle, locVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, *locCovarianceMatrix))
		return false;

	locKinematicData->setMomentum(DVector3(locMomentum.X(),locMomentum.Y(),locMomentum.Z()));
	locKinematicData->setPosition(DVector3(locSpacetimeVertex.Vect().X(),locSpacetimeVertex.Vect().Y(),locSpacetimeVertex.Vect().Z()));
	locKinematicData->setTime(locSpacetimeVertex.T());
	locKinematicData->setErrorMatrix(*locCovarianceMatrix);
	locKinematicData->setPathLength(locPathLengthPair.first, locPathLengthPair.second);
	return true;
}
