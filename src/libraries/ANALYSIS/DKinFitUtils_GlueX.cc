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
}

void DKinFitUtils_GlueX::Set_MaxPoolSizes(size_t locNumReactions, size_t locExpectedNumCombos)
{
	Set_MaxKinFitConstraintVertexPoolSize(2*2*locNumReactions*locExpectedNumCombos*2); //extra x2: guess fit!
	Set_MaxKinFitConstraintSpacetimePoolSize(2*2*locNumReactions*locExpectedNumCombos*2); //extra x2: guess fit!
	Set_MaxKinFitConstraintP4PoolSize(locNumReactions*locExpectedNumCombos*2);
	Set_MaxKinFitConstraintMassPoolSize(2*locNumReactions*locExpectedNumCombos*2);

	Set_MaxKinFitChainPoolSize(locNumReactions*locExpectedNumCombos);
	Set_MaxKinFitChainStepPoolSize(3*locNumReactions*locExpectedNumCombos);
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

	DKinFitUtils::Reset_NewEvent();
}

bool DKinFitUtils_GlueX::Get_IncludeBeamlineInVertexFitFlag(void) const
{
	return dIncludeBeamlineInVertexFitFlag; //at least until covariance matrix is set for beam photons
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

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_BeamParticle(const DBeamPhoton* locBeamPhoton)
{
	pair<const DBeamPhoton*, const DEventRFBunch*> locSourcePair(locBeamPhoton, NULL);
	if(dParticleMap_SourceToInput_Beam.find(locSourcePair) != dParticleMap_SourceToInput_Beam.end())
		return dParticleMap_SourceToInput_Beam[locSourcePair]; //not unique, return existing

	TLorentzVector locSpacetimeVertex(Make_TVector3(locBeamPhoton->position()), locBeamPhoton->time());
	TVector3 locMomentum = Make_TVector3(locBeamPhoton->momentum());
	Particle_t locPID = locBeamPhoton->PID();

	auto locKinFitParticle = DKinFitUtils::Make_BeamParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, locBeamPhoton->errorMatrix());
	dParticleMap_SourceToInput_Beam[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locBeamPhoton;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_BeamParticle(const DBeamPhoton* locBeamPhoton, const DEventRFBunch* locEventRFBunch)
{
	pair<const DBeamPhoton*, const DEventRFBunch*> locSourcePair(locBeamPhoton, locEventRFBunch);
	if(dParticleMap_SourceToInput_Beam.find(locSourcePair) != dParticleMap_SourceToInput_Beam.end())
		return dParticleMap_SourceToInput_Beam[locSourcePair]; //not unique, return existing

	//set rf time for beam particle
	TLorentzVector locSpacetimeVertex(Make_TVector3(locBeamPhoton->position()), locEventRFBunch->dTime);
	TVector3 locMomentum = Make_TVector3(locBeamPhoton->momentum());
	Particle_t locPID = locBeamPhoton->PID();

	//set rf time variance in covariance matrix
	auto locCovarianceMatrix = Get_SymMatrixResource(7);
	*locCovarianceMatrix = *(locBeamPhoton->errorMatrix());
	(*locCovarianceMatrix)(6, 6) = locEventRFBunch->dTimeVariance;
	//zero the correlation terms
	for(int loc_i = 0; loc_i < 6; ++loc_i)
	{
		(*locCovarianceMatrix)(6, loc_i) = 0.0;
		(*locCovarianceMatrix)(loc_i, 6) = 0.0;
	}

	auto locKinFitParticle = DKinFitUtils::Make_BeamParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, locCovarianceMatrix);
	dParticleMap_SourceToInput_Beam[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locBeamPhoton;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_DetectedParticle(const DKinematicData* locKinematicData)
{
	if(dParticleMap_SourceToInput_DetectedParticle.find(locKinematicData) != dParticleMap_SourceToInput_DetectedParticle.end())
		return dParticleMap_SourceToInput_DetectedParticle[locKinematicData]; //not unique, return existing

	TLorentzVector locSpacetimeVertex(Make_TVector3(locKinematicData->position()), locKinematicData->time());
	TVector3 locMomentum = Make_TVector3(locKinematicData->momentum());
	Particle_t locPID = locKinematicData->PID();

	double locPathLength = 0.0;
	auto locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locKinematicData);
	if(locChargedHypo != nullptr)
		locPathLength = locChargedHypo->Get_PathLength();
	else
	{
		auto locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locKinematicData);
		if(locNeutralHypo != nullptr)
			locPathLength = locNeutralHypo->Get_PathLength();
	}

	auto locKinFitParticle = DKinFitUtils::Make_DetectedParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, locPathLength, locKinematicData->errorMatrix());
	dParticleMap_SourceToInput_DetectedParticle[locKinematicData] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locKinematicData;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_DetectedShower(const DNeutralShower* locNeutralShower, Particle_t locPID)
{
	pair<const DNeutralShower*, Particle_t> locSourcePair(locNeutralShower, locPID);
	if(dParticleMap_SourceToInput_Shower.find(locSourcePair) != dParticleMap_SourceToInput_Shower.end())
		return dParticleMap_SourceToInput_Shower[locSourcePair]; //not unique, return existing

	//use DNeutralShower object (doesn't make assumption about vertex!)
	TLorentzVector locShowerSpacetime = Make_TLorentzVector(locNeutralShower->dSpacetimeVertex);
	auto locKinFitParticle = DKinFitUtils::Make_DetectedShower(PDGtype(locPID), ParticleMass(locPID), locShowerSpacetime, locNeutralShower->dEnergy, locNeutralShower->dCovarianceMatrix);

	dParticleMap_SourceToInput_Shower[locSourcePair] = locKinFitParticle;
	dParticleMap_InputToSource_JObject[locKinFitParticle] = locNeutralShower;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_TargetParticle(Particle_t locPID)
{
	if(dParticleMap_SourceToInput_Target.find(locPID) != dParticleMap_SourceToInput_Target.end())
		return dParticleMap_SourceToInput_Target[locPID]; //not unique, return existing

	auto locKinFitParticle = DKinFitUtils::Make_TargetParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMap_SourceToInput_Target[locPID] = locKinFitParticle;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_MissingParticle(Particle_t locPID)
{
	if(dParticleMap_SourceToInput_Missing.find(locPID) != dParticleMap_SourceToInput_Missing.end())
		return dParticleMap_SourceToInput_Missing[locPID]; //not unique, return existing

	auto locKinFitParticle = DKinFitUtils::Make_MissingParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMap_SourceToInput_Missing[locPID] = locKinFitParticle;
	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils_GlueX::Make_DecayingParticle(Particle_t locPID, const set<shared_ptr<DKinFitParticle>>& locFromInitialState, const set<shared_ptr<DKinFitParticle>>& locFromFinalState)
{
	DDecayingParticleInfo locDecayingParticleInfo(locPID, locFromInitialState, locFromFinalState);
	if(dParticleMap_SourceToInput_Decaying.find(locDecayingParticleInfo) != dParticleMap_SourceToInput_Decaying.end())
		return dParticleMap_SourceToInput_Decaying[locDecayingParticleInfo]; //not unique, return existing

	auto locKinFitParticle = DKinFitUtils::Make_DecayingParticle(PDGtype(locPID), ParticleCharge(locPID), ParticleMass(locPID), locFromInitialState, locFromFinalState);
	dParticleMap_SourceToInput_Decaying[locDecayingParticleInfo] = locKinFitParticle;
	dParticleMap_InputToSource_Decaying.emplace(locKinFitParticle, locDecayingParticleInfo);
	return locKinFitParticle;
}

/**************************************************************** MAKE DKINFITCHAIN ****************************************************************/

const DKinFitChain* DKinFitUtils_GlueX::Make_KinFitChain(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType)
{
	//locKinFitType input in case want to do a different fit
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: Create DKinFitChain." << endl;

	DKinFitChain* locKinFitChain = Get_KinFitChainResource();
	locKinFitChain->Set_DefinedParticleStepIndex(-1); //unless changed below

	//Make chain, excluding decaying particles
		//They must be created using the detected particles, so just create those first
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		//Create it
		auto locKinFitChainStep = Make_KinFitChainStep(locReactionVertexInfo, locReaction, locParticleCombo, locKinFitType, loc_i, locKinFitChain);
		locKinFitChain->Add_KinFitChainStep(locKinFitChainStep);
	}

	//Now define the decaying particles using the detected & beam particles

	//start looping from the back, using invariant mass wherever possible
		//not possible if missing decay product: will then compute via missing mass in the next step
	//but first, figure out which steps we must use missing mass for (and thus must skip on this pass below)
	set<size_t> locMissingMassSteps;
	auto locDefinedParticleStepIndices = DAnalysis::Get_DefinedParticleStepIndex(locReaction);
	for(int locCurrentStepIndex : locDefinedParticleStepIndices)
	{
		while(locCurrentStepIndex != -1)
		{
			if((locCurrentStepIndex == 0) && locKinFitChain->Get_KinFitChainStep(0)->Get_InitialParticles().empty())
				break; //is an open-ended decaying particle in the initial state: is ok to define via invariant mass
			locMissingMassSteps.insert(locCurrentStepIndex);
			locCurrentStepIndex = locKinFitChain->Get_KinFitChainStep(locCurrentStepIndex)->Get_InitialParticleDecayFromStepIndex();
		}
	}

	//now create decaying particles
	//do the invariant mass loop
	set<size_t> locProcessedSteps;
	for(int loc_i = locParticleCombo->Get_NumParticleComboSteps() - 1; loc_i >= 0; --loc_i)
	{
		//get steps
		auto locReactionStep = locReaction->Get_ReactionStep(loc_i);
		auto locKinFitChainStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(loc_i));

		//skip steps that must defined via missing mass
		if(locMissingMassSteps.find(loc_i) != locMissingMassSteps.end())
			continue; //decay products contain a missing or open-ended-decaying particle

		//check if initial particle already created for this step (i.e. beam)
		if(locKinFitChainStep->Get_InitialParticle(0) != nullptr)
			continue; //only do decaying particles

		//mark the progress
		locProcessedSteps.insert(loc_i);

		//don't create particle for omega, etc. (will confuse kinfitter)
		auto locPID = locReactionStep->Get_InitialPID();
		if(!IsFixedMass(locPID))
			continue;

		//since going in reverse order, all decay products are ready: create the decaying particle
		auto locDecaySourceParticles = Get_StepParticles_NonNull(locKinFitChain, locReaction, loc_i, DReactionStep::Get_ParticleIndex_Initial());
		auto locDecayingParticle = Make_DecayingParticle(locPID, locDecaySourceParticles.first, locDecaySourceParticles.second);

		//set decaying particle in the chain
		locKinFitChainStep->Set_InitialParticle(locDecayingParticle, 0);
		auto locFromIndices = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, loc_i);
		if((locFromIndices.first < 0) || (locFromIndices.first == loc_i))
			continue; //intial-state open-ended decaying particle

		auto locProductionStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locFromIndices.first));
		locProductionStep->Set_FinalParticle(locDecayingParticle, locFromIndices.second);
		locKinFitChain->Set_DecayStepIndex(locDecayingParticle, loc_i);
	}

	//now loop from the front, using missing mass for the rest
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		auto locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locProcessedSteps.find(loc_i) != locProcessedSteps.end())
			continue;

		auto locReactionStep = locReaction->Get_ReactionStep(loc_i);
		auto locKinFitChainStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(loc_i));

		//create decaying particle in final state: first figure out which one needs a particle
		for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
		{
			//DecayStepIndex: >= 0 if decaying, where the # is the step representing the particle decay
			int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, loc_i, loc_j);
			if(locDecayStepIndex < 0)
				continue; //not decaying

			if(locKinFitChain->Get_KinFitChainStep(locDecayStepIndex)->Get_InitialParticle(0) != nullptr)
				continue; //decaying particle already created

			//don't create particle for omega, etc. (will confuse kinfitter)
			auto locPID = locReactionStep->Get_FinalPID(loc_j);
			if(!IsFixedMass(locPID))
				continue;

			//make decaying particle
			auto locDecaySourceParticles = Get_StepParticles_NonNull(locKinFitChain, locReaction, loc_i, loc_j);
			auto locDecayingParticle = Make_DecayingParticle(locPID, locDecaySourceParticles.first, locDecaySourceParticles.second);
			locKinFitChainStep->Set_FinalParticle(locDecayingParticle, loc_j);

			auto locDecayStep = const_cast<DKinFitChainStep*>(locKinFitChain->Get_KinFitChainStep(locDecayStepIndex));
			locDecayStep->Set_InitialParticle(locDecayingParticle, 0);
			locKinFitChain->Set_DecayStepIndex(locDecayingParticle, locDecayStepIndex);
		}
	}

	if(dDebugLevel > 10)
	{
		cout << "DKinFitUtils_GlueX: DKinFitChain Created. Printing:" << endl;
		locKinFitChain->Print_InfoToScreen();
	}

	return locKinFitChain;
}

pair<set<shared_ptr<DKinFitParticle>>, set<shared_ptr<DKinFitParticle>>> DKinFitUtils_GlueX::Get_StepParticles_NonNull(const DKinFitChain* locKinFitChain, const DReaction* locReaction, size_t locStepIndex, int locNonFixedMassParticleIndex) const
{
	//If one of the final particles is null (decaying non-fixed-mass particle), remove it, AND go to its decay step and get ALL particles (put in init or final state)
	auto locKinFitStep = locKinFitChain->Get_KinFitChainStep(locStepIndex);

	//Get Initial particles, erase nulls
	auto locGrabbedInitialParticles = locKinFitStep->Get_InitialParticles();
	set<shared_ptr<DKinFitParticle>> locInitialParticles(locGrabbedInitialParticles.begin(), locGrabbedInitialParticles.end());
	locInitialParticles.erase(nullptr); //remove any null particles

	//Get Final particles, erase nulls
	auto locGrabbedFinalParticles = locKinFitStep->Get_FinalParticles();
	set<shared_ptr<DKinFitParticle>> locFinalParticles(locGrabbedFinalParticles.begin(), locGrabbedFinalParticles.end());
	locFinalParticles.erase(nullptr); //remove any null particles

	//If nulls in initial state: go to the previous step
	for(size_t loc_i = 0; loc_i < locGrabbedInitialParticles.size(); ++loc_i)
	{
		if((locGrabbedInitialParticles[loc_i] != nullptr) || (locNonFixedMassParticleIndex == DReactionStep::Get_ParticleIndex_Initial()))
			continue;

		//null particles in initial state: go to the previous step
		auto locParticlePair = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locStepIndex);
		auto locStepParticles = Get_StepParticles_NonNull(locKinFitChain, locReaction, locParticlePair.first, locParticlePair.second);
		locInitialParticles.insert(locStepParticles.first.begin(), locStepParticles.first.end());
		locFinalParticles.insert(locStepParticles.second.begin(), locStepParticles.second.end());
	}

	//If nulls in final state: go to the next step
	for(size_t loc_i = 0; loc_i < locGrabbedFinalParticles.size(); ++loc_i)
	{
		if((locGrabbedFinalParticles[loc_i] != nullptr) || (locNonFixedMassParticleIndex == int(loc_i)))
			continue;

		auto locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		auto locStepParticles = Get_StepParticles_NonNull(locKinFitChain, locReaction, locDecayStepIndex, DReactionStep::Get_ParticleIndex_Initial());
		locInitialParticles.insert(locStepParticles.first.begin(), locStepParticles.first.end());
		locFinalParticles.insert(locStepParticles.second.begin(), locStepParticles.second.end());
	}

	if(dDebugLevel >= 10)
	{
		cout << "DKinFitUtils_GlueX::Get_StepParticles_NonNull:" << endl;
		cout << "reaction, step-index, non-fixed-mass particle index: " << locReaction->Get_ReactionName() << ", " << locStepIndex << ", " << locNonFixedMassParticleIndex << endl;
		cout << "Chain: " << endl;
		locKinFitChain->Print_InfoToScreen();
		cout << "Init-particles:" << endl;
		for(auto& locParticle : locInitialParticles)
			cout << locParticle->Get_PID() << ", " << locParticle << endl;
		cout << "Final-particles:" << endl;
		for(auto& locParticle : locFinalParticles)
			cout << locParticle->Get_PID() << ", " << locParticle << endl;
	}
	return std::make_pair(locInitialParticles, locFinalParticles);
}

DKinFitChainStep* DKinFitUtils_GlueX::Make_KinFitChainStep(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, DKinFitType locKinFitType, size_t locStepIndex, DKinFitChain* locKinFitChain)
{
	//Start and register new step
	auto locKinFitChainStep = Get_KinFitChainStepResource();
	locKinFitChainStep->Set_ConstrainDecayingMassFlag(false); //unless changed later

	//get the steps
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(locStepIndex);
	const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	auto locStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(locStepIndex);
	int locKinFitStepIndex = locKinFitChain->Get_NumKinFitChainSteps();

	//if doing a vertex fit, see which neutral particles can be treated as showers
	bool locSpactimeIsFitFlag = (locKinFitType == d_SpacetimeFit) || (locKinFitType == d_P4AndSpacetimeFit);

	//initial beam particle
	locKinFitChainStep->Set_InitialParticleDecayFromStepIndex(-1); //unless set otherwise below (enclosed decaying particle)
	const DBeamPhoton* locBeamPhoton = dynamic_cast<const DBeamPhoton*>(locParticleComboStep->Get_InitialParticle_Measured());
	if(locBeamPhoton != NULL)
		locKinFitChainStep->Add_InitialParticle(Make_BeamParticle(locBeamPhoton));
	else //decaying particle
	{
		locKinFitChainStep->Add_InitialParticle(nullptr); //will change later
		if(locStepIndex == 0) //open-ended
			locKinFitChain->Set_DefinedParticleStepIndex(locKinFitStepIndex);
		else //enclosed
		{
			int locDecayFromStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locStepIndex).first;
			locKinFitChainStep->Set_InitialParticleDecayFromStepIndex(locDecayFromStepIndex);
			Particle_t locPID = locReactionStep->Get_InitialPID();
			auto locConstrainMassFlag = IsFixedMass(locPID) ? locReactionStep->Get_KinFitConstrainInitMassFlag() : false;
			locKinFitChainStep->Set_ConstrainDecayingMassFlag(locConstrainMassFlag);
		}
		//will create decaying particle later
	}

	//target particle
	Particle_t locTargetPID = locReactionStep->Get_TargetPID();
	if(locTargetPID != Unknown)
		locKinFitChainStep->Add_InitialParticle(Make_TargetParticle(locTargetPID));

	//final state particles
	for(size_t loc_j = 0; loc_j < locParticleComboStep->Get_NumFinalParticles(); ++loc_j)
	{
		int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_j);
		const DKinematicData* locKinematicData = locParticleComboStep->Get_FinalParticle_Measured(loc_j);
		Particle_t locPID = locReactionStep->Get_FinalPID(loc_j);

		if(locReactionStep->Get_MissingParticleIndex() == int(loc_j)) //missing particle
		{
			locKinFitChain->Set_DefinedParticleStepIndex(locKinFitStepIndex);
			locKinFitChainStep->Add_FinalParticle(Make_MissingParticle(locPID));
		}
		else if(locDecayStepIndex >= 0) //decaying particle
			locKinFitChainStep->Add_FinalParticle(nullptr); //skip for now, will create later (unless not fixed mass)
		else if(ParticleCharge(locPID) == 0) //detected neutral
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locKinematicData);

			//Determine whether we should use the particle or the shower object
			bool locNeutralShowerFlag = locStepVertexInfo->Get_FittableVertexFlag();
			if((ParticleMass(locPID) > 0.0) && !locSpactimeIsFitFlag)
				locNeutralShowerFlag = false; //massive shower momentum is defined by t, which isn't fit: use particle

			if(!locNeutralShowerFlag)
				locKinFitChainStep->Add_FinalParticle(Make_DetectedParticle(locNeutralParticleHypothesis));
			else //in a vertex constraint: make shower
			{
				const DNeutralShower* locNeutralShower = locNeutralParticleHypothesis->Get_NeutralShower();
				locKinFitChainStep->Add_FinalParticle(Make_DetectedShower(locNeutralShower, locNeutralParticleHypothesis->PID()));
			}
		}
		else //detected charged track
		{
			const DChargedTrackHypothesis* locChargedTrackHypothesis = static_cast<const DChargedTrackHypothesis*>(locKinematicData);
			locKinFitChainStep->Add_FinalParticle(Make_DetectedParticle(locChargedTrackHypothesis));
		}
	}

	//Inclusive?
	if(locReactionStep->Get_MissingParticleIndex() == DReactionStep::Get_ParticleIndex_Inclusive())
		locKinFitChain->Set_IsInclusiveChannelFlag(true);

	return locKinFitChainStep;
}

/**************************************************************** MAKE CONSTRAINTS *****************************************************************/

set<DKinFitConstraint*> DKinFitUtils_GlueX::Create_Constraints(const DReactionVertexInfo* locReactionVertexInfo, const DReaction* locReaction, const DParticleCombo* locParticleCombo, const DKinFitChain* locKinFitChain, DKinFitType locKinFitType, deque<DKinFitConstraint_Vertex*>& locSortedVertexConstraints)
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: Create constraints." << endl;

	//All constraints
	set<DKinFitConstraint*> locAllConstraints;

	//Create Mass Constraints
	set<DKinFitConstraint_Mass*> locMassConstraints;
	map<shared_ptr<DKinFitParticle>, DKinFitConstraint_Mass*> locParticleMassConstraintMap;
	map<shared_ptr<DKinFitParticle>, size_t> locParticleDecayStepMap; //key is decaying particle, value is step index
	map<size_t, DKinFitConstraint_Mass*> locStepMassConstraintMap;
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit))
	{
		//loop over steps, but skip init step (if open-ended decaying, do p4 constraint instead)
		for(size_t loc_i = 1; loc_i < locKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
		{
			auto locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(loc_i);
			if(!locKinFitChainStep->Get_ConstrainDecayingMassFlag())
				continue; //don't apply mass constraint to this step

			auto locInitialParticles = locKinFitChainStep->Get_InitialParticles();
			auto locParticleIterator = locInitialParticles.begin();
			for(; locParticleIterator != locInitialParticles.end(); ++locParticleIterator)
			{
				if((*locParticleIterator) == nullptr)
					continue;
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
		//don't do if inclusive reaction, or more than one missing particle
		//pick the step containing the defined (missing or open-ended-decaying) particle
		//if no defined particle, use the first step
	DKinFitConstraint_P4* locP4Constraint = NULL;
	auto locDefinedParticleStepIndices = DAnalysis::Get_DefinedParticleStepIndex(locReaction);
	if(!locKinFitChain->Get_IsInclusiveChannelFlag() && (locDefinedParticleStepIndices.size() <= 1) && ((locKinFitType == d_P4Fit) || (locKinFitType == d_P4AndVertexFit) || (locKinFitType == d_P4AndSpacetimeFit)))
	{
		int locDefinedParticleStepIndex = locDefinedParticleStepIndices.empty() ? -1 : locDefinedParticleStepIndices[0];
		int locP4StepIndex = (locDefinedParticleStepIndex >= 0) ? locDefinedParticleStepIndex : 0;
		auto locStepParticles = Get_StepParticles_NonNull(locKinFitChain, locReaction, locP4StepIndex);
		locP4Constraint = Make_P4Constraint(locStepParticles.first, locStepParticles.second);

		//OK, now, check to see if the system is overly constrained: 
			//there must be at least one particle with non-zero errors in the p4 constraint that is NOT in a mass constraint
		bool locNonZeroErrorFlag = false;
		set<size_t> locP4ConstrainedParticleSteps;

		auto locAllParticles = locP4Constraint->Get_AllParticles();
		for(auto& locKinFitParticle : locAllParticles)
		{
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
					locDefinedP4 = dAnalysisUtilities->Calc_MissingP4(locReaction, locParticleCombo, false);
				else
					locDefinedP4 = dAnalysisUtilities->Calc_FinalStateP4(locReaction, locParticleCombo, 0, false);
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
		for(auto& locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		{
			if(!locStepVertexInfo->Get_FittableVertexFlag())
				continue;
			auto locDX4 = locParticleCombo->Get_ParticleComboStep(locStepVertexInfo->Get_StepIndices().front())->Get_SpacetimeVertex();
			TLorentzVector locX4(locDX4.X(), locDX4.Y(), locDX4.Z(), locDX4.T());

			auto locFullConstrainParticles = Build_ParticleSet(locStepVertexInfo->Get_FullConstrainParticles(true), locKinFitChain);
			auto locNoConstrainParticles = Build_ParticleSet(locStepVertexInfo->Get_NoConstrainParticles(true), locKinFitChain);
			auto locOnlyConstrainTimeParticles = Build_ParticleSet(locStepVertexInfo->Get_OnlyConstrainTimeParticles(), locKinFitChain);
			if(locSpacetimeFitFlag)
				locSortedVertexConstraints.push_back(Make_SpacetimeConstraint(locFullConstrainParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles, locX4));
			else
			{
				locNoConstrainParticles.insert(locOnlyConstrainTimeParticles.begin(), locOnlyConstrainTimeParticles.end());
				locSortedVertexConstraints.push_back(Make_VertexConstraint(locFullConstrainParticles, locNoConstrainParticles, locX4.Vect()));
			}
		}
	}
	locAllConstraints.insert(locSortedVertexConstraints.begin(), locSortedVertexConstraints.end());

	if(dDebugLevel > 10)
		cout << "DKinFitUtils_GlueX: All Constraints Created." << endl;

	return locAllConstraints;
}

set<shared_ptr<DKinFitParticle>> DKinFitUtils_GlueX::Build_ParticleSet(const vector<pair<int, int>>& locParticleIndices, const DKinFitChain* locKinFitChain)
{
	set<shared_ptr<DKinFitParticle>> locParticles;
	for(auto& locParticlePair : locParticleIndices)
	{
		auto locStep = locKinFitChain->Get_KinFitChainStep(locParticlePair.first);
		if(locParticlePair.second >= 0)
			locParticles.insert(locStep->Get_FinalParticle(locParticlePair.second));
		else if(locParticlePair.second == DReactionStep::Get_ParticleIndex_Initial())
			locParticles.insert(locStep->Get_InitialParticle(0));
		else
			locParticles.insert(locStep->Get_InitialParticle(1));
	}
	return locParticles;
}

/************************************************************** CONSTRAINT PREDICTORS **************************************************************/

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

			Particle_t locPID = locReactionStep->Get_InitialPID();
			if(!IsFixedMass(locPID))
				continue; //don't apply mass constraint to this step

			locMassConstraintStrings[loc_i] = string("#it{m}_{") + string(ParticleName_ROOT(locPID)) + string("}");
		}

		//P4 Constraint //don't do if inclusive reaction, or more than 1 missing particle
		auto locDefinedParticleStepIndices = DAnalysis::Get_DefinedParticleStepIndex(locReaction);
		if(!locReaction->Get_IsInclusiveFlag() && (locDefinedParticleStepIndices.size() <= 1))
		{
			int locDefinedParticleStepIndex = locDefinedParticleStepIndices.empty() ? -1 : locDefinedParticleStepIndices[0];
			//OK, now, check to see if the system is overly constrained: 
				//there must be at least one particle with non-zero errors in the p4 constraint that is NOT in a mass constraint
			bool locNonZeroErrorFlag = false;

			//Find step used for p4 constraint: the one with missing/open-ended-decaying particle (if any: else first step)
			int locP4StepIndex = (locDefinedParticleStepIndex >= 0) ? locDefinedParticleStepIndex : 0;
			const DReactionStep* locReactionStep = locReaction->Get_ReactionStep(locP4StepIndex);

			set<size_t> locP4ConstrainedParticleSteps;

			//check if initial & final states if have non-zero cov errors:
			if((locP4StepIndex == 0) && (locReactionStep->Get_TargetPID() != Unknown) && dWillBeamHaveErrorsFlag)
				locNonZeroErrorFlag = true; //beam: we're good
			else if((locP4StepIndex != 0) && (locMassConstraintStrings.find(locP4StepIndex) == locMassConstraintStrings.end()))
				locNonZeroErrorFlag = true; //decaying particle, but mass not constrained: we're good (unless it's e.g. an omega. ugh.)
			else //check final state
			{
				//check if final state has non-zero cov errors (detected):
				for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalPIDs(); ++loc_i)
				{
					if(locReactionStep->Get_MissingParticleIndex() == int(loc_i))
						continue; //missing

					int locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locP4StepIndex, loc_i);
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
		auto locConstraintTuple = Predict_VertexConstraints(locReactionVertexInfo, locSpacetimeFitFlag);
		if(std::get<0>(locConstraintTuple) > 0)
		{
			locNumConstraints += std::get<0>(locConstraintTuple);
			locNumUnknowns += std::get<1>(locConstraintTuple);
			if(locAllConstraintsString != "")
				locAllConstraintsString += ", ";
			locAllConstraintsString += std::get<2>(locConstraintTuple);
		}
	}

	return locAllConstraintsString;
}

tuple<size_t, size_t, string> DKinFitUtils_GlueX::Predict_VertexConstraints(const DReactionVertexInfo* locReactionVertexInfo, bool locSpacetimeFitFlag) const
{
	//returned: #constraints, constraint string
	size_t locNumConstraints = 0, locNumUnknowns = 0;
	string locAllConstraintString;
	auto locStepVertexInfos = locReactionVertexInfo->Get_StepVertexInfos();
	for(auto& locVertexInfo : locStepVertexInfos)
	{
		if(!locVertexInfo->Get_FittableVertexFlag())
			continue;

		auto locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles(true);
		if(locSpacetimeFitFlag)
		{
			locNumConstraints += 3*locFullConstrainParticles.size();
			locNumConstraints += locVertexInfo->Get_OnlyConstrainTimeParticles().size();
			locNumUnknowns += 4;
		}
		else //vertex only
		{
			locNumConstraints += 2*locFullConstrainParticles.size();
			locNumUnknowns += 3;
		}

		//add to the full constraint string
		if(locAllConstraintString != "")
			locAllConstraintString += ", ";
		locAllConstraintString += Build_VertexConstraintString(locVertexInfo, locSpacetimeFitFlag);
	}

	return std::make_tuple(locNumConstraints, locNumUnknowns, locAllConstraintString);
}

string DKinFitUtils_GlueX::Build_VertexConstraintString(const DReactionStepVertexInfo* locVertexInfo, bool locSpacetimeFitFlag) const
{
	auto locReaction = locVertexInfo->Get_Reaction();
	string locConstraintString = locSpacetimeFitFlag ? "#it{x}^{4}_{" : "#it{x}^{3}_{";

	//if a decaying particle decays in-place at this vertex, we don't want to put it into the string
	//this will happen if: the vertex-info represents more than one step
	//they will be present in the info as non-constrain particles in the final-state
	auto locStepIndices = locVertexInfo->Get_StepIndices(); //steps are listed in order
	set<pair<int, int>> locExcludeDecayParticleIndices; //these are final-state indices
	for(size_t loc_i = 1; loc_i < locStepIndices.size(); ++loc_i)
		locExcludeDecayParticleIndices.insert(DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locStepIndices[loc_i]));

	//initial state
	auto locParticles = locVertexInfo->Get_Particles(d_InitialState);
	auto locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles(true, d_InitialState);
	for(auto locIndices : locParticles)
	{
		auto locStep = locReaction->Get_ReactionStep(locIndices.first);
		Particle_t locPID = locStep->Get_PID(locIndices.second);
		string locParticleString = ParticleName_ROOT(locPID);
		if(locIndices.second == locStep->Get_MissingParticleIndex())
			locConstraintString += string("(") + locParticleString + string(")"); //missing
		else if(std::binary_search(locFullConstrainParticles.begin(), locFullConstrainParticles.end(), locIndices)) //constraining
			locConstraintString += string("#color[4]{") + locParticleString + string("}"); //blue
		else //no-constrain
			locConstraintString += locParticleString; //plain
	}

	//final state
	locConstraintString += "#rightarrow";
	locParticles = locVertexInfo->Get_Particles(d_FinalState);
	locFullConstrainParticles = locVertexInfo->Get_FullConstrainParticles(true, d_FinalState);
	auto locOnlyConstrainTimeParticles = locVertexInfo->Get_OnlyConstrainTimeParticles();
	for(auto locIndices : locParticles)
	{
		if(locExcludeDecayParticleIndices.find(locIndices) != locExcludeDecayParticleIndices.end())
			continue; //exclude no-constrain decaying particle

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
	auto locCovarianceMatrix = Get_SymMatrixResource(7);
	if(!DKinFitUtils::Propagate_TrackInfoToCommonVertex(locKinFitParticle, locVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, locCovarianceMatrix.get()))
		return false;

	locKinematicData->setMomentum(DVector3(locMomentum.X(),locMomentum.Y(),locMomentum.Z()));
	locKinematicData->setPosition(DVector3(locSpacetimeVertex.Vect().X(),locSpacetimeVertex.Vect().Y(),locSpacetimeVertex.Vect().Z()));
	locKinematicData->setTime(locSpacetimeVertex.T());
	locKinematicData->setErrorMatrix(locCovarianceMatrix);
	return true;
}
