#include "ANALYSIS/DParticleComboCreator.h"
#include "ANALYSIS/DSourceComboer.h"
#include "ANALYSIS/DAnalysisUtilities.h"

namespace DAnalysis
{

DParticleComboCreator::DParticleComboCreator(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboTimeHandler* locSourceComboTimeHandler, const DSourceComboVertexer* locSourceComboVertexer) :
		dSourceComboer(locSourceComboer), dSourceComboTimeHandler(locSourceComboTimeHandler), dSourceComboVertexer(locSourceComboVertexer)
{
	gPARMS->SetDefaultParameter("COMBO:DEBUG_LEVEL", dDebugLevel);
	dKinFitUtils = new DKinFitUtils_GlueX(locEventLoop);

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses); //make sure that brun() is called for the default factory!!!
	dNeutralParticleHypothesisFactory = static_cast<DNeutralParticleHypothesis_factory*>(locEventLoop->GetFactory("DNeutralParticleHypothesis"));

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses); //make sure that brun() is called for the default factory!!!
	dChargedTrackHypothesisFactory = static_cast<DChargedTrackHypothesis_factory*>(locEventLoop->GetFactory("DChargedTrackHypothesis"));

	vector<const DBeamPhoton*> locBeamPhotons;
	locEventLoop->Get(locBeamPhotons); //make sure that brun() is called for the default factory!!!
	dBeamPhotonfactory = static_cast<DBeamPhoton_factory*>(locEventLoop->GetFactory("DBeamPhoton"));

	locEventLoop->GetSingle(dParticleID);

	//error matrix //too lazy to compute properly right now ... need to hack DAnalysisUtilities::Calc_DOCA()
	dVertexCovMatrix.ResizeTo(4, 4);
	dVertexCovMatrix.Zero();
	dVertexCovMatrix(0, 0) = 0.65; //x variance //from monitoring plots of vertex
	dVertexCovMatrix(1, 1) = 0.65; //y variance //from monitoring plots of vertex
	dVertexCovMatrix(2, 2) = 1.5; //z variance //a guess, semi-guarding against the worst case scenario //ugh
	dVertexCovMatrix(3, 3) = 0.0; //t variance //not used
}

void DParticleComboCreator::Reset(void)
{
	for(const auto& locRFPair : dRFBunchMap)
		dResourcePool_EventRFBunch.Recycle(locRFPair.second);
	dRFBunchMap.clear();

	for(const auto& locStepPair : dComboStepMap)
		dResourcePool_ParticleComboStep.Recycle(locStepPair.second);
	dComboStepMap.clear();

	for(const auto& locComboPair : dComboMap)
		dResourcePool_ParticleCombo.Recycle(locComboPair.second);
	dComboMap.clear();

	for(const auto& locHypoPair : dChargedHypoMap)
		dChargedTrackHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dChargedHypoMap.clear();

	for(const auto& locHypoPair : dKinFitChargedHypoMap)
		dChargedTrackHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dKinFitChargedHypoMap.clear();

	for(const auto& locHypoPair : dNeutralHypoMap)
		dNeutralParticleHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dNeutralHypoMap.clear();

	for(const auto& locHypoPair : dKinFitNeutralHypoMap)
		dNeutralParticleHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dKinFitNeutralHypoMap.clear();

	for(const auto& locBeamPair : dKinFitBeamPhotonMap)
		dBeamPhotonfactory->Recycle_Resource(locBeamPair.second);
	dKinFitBeamPhotonMap.clear();
}

bool DParticleComboCreator::Get_CreateNeutralErrorMatrixFlag_Combo(const DReactionVertexInfo* locReactionVertexInfo, DKinFitType locKinFitType)
{
	//is there at least one non-fittable vertex that has neutrals?
	auto locDanglingIterator = dDanglingNeutralsFlagMap.find(locReactionVertexInfo);
	bool locDanglingNeutralsFlag = false;
	if(locDanglingIterator != dDanglingNeutralsFlagMap.end())
		locDanglingNeutralsFlag = locDanglingIterator->second;
	else
	{
		for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		{
			if(locStepVertexInfo->Get_FittableVertexFlag())
				continue;
			if(!locStepVertexInfo->Get_OnlyConstrainTimeParticles().empty())
			{
				locDanglingNeutralsFlag = true; //photons
				break;
			}
			if(!locStepVertexInfo->Get_NoConstrainParticles(true, d_FinalState, d_Neutral, false, false, false).empty())
			{
				locDanglingNeutralsFlag = true; //massive neutrals
				break;
			}
		}
		dDanglingNeutralsFlagMap.emplace(locReactionVertexInfo, locDanglingNeutralsFlag);
	}
	return ((locKinFitType != d_NoFit) && ((locKinFitType == d_P4Fit) || locDanglingNeutralsFlag));
}

const DParticleCombo* DParticleComboCreator::Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift, DKinFitType locKinFitType)
{
	if(dDebugLevel > 0)
		cout << "Building particle combo" << endl;
	auto locCreateNeutralErrorMatrixFlag_Combo = Get_CreateNeutralErrorMatrixFlag_Combo(locReactionVertexInfo, locKinFitType);
	auto locComboTuple = std::make_tuple(locReactionVertexInfo, locFullCombo, locBeamParticle, locRFBunchShift, locCreateNeutralErrorMatrixFlag_Combo);
	auto locComboIterator = dComboMap.find(locComboTuple);
	if(locComboIterator != dComboMap.end())
	{
		if(dDebugLevel > 0)
			cout << "Combo previously created, returning." << endl;
		return locComboIterator->second;
	}

	auto locParticleCombo = dResourcePool_ParticleCombo.Get_Resource();
	locParticleCombo->Reset();
	dComboMap.emplace(locComboTuple, locParticleCombo);

	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locFullCombo, locBeamParticle).Z();

	//Get/Create RF Bunch
	const DEventRFBunch* locEventRFBunch = nullptr;
	auto locRFIterator = dRFBunchMap.find(locRFBunchShift);
	if(locRFIterator != dRFBunchMap.end())
		locEventRFBunch = locRFIterator->second;
	else
	{
		auto locInitialEventRFBunch = dSourceComboTimeHandler->Get_InitialEventRFBunch();
		auto locNewEventRFBunch = dResourcePool_EventRFBunch.Get_Resource();
		auto locRFTime = dSourceComboTimeHandler->Calc_RFTime(locRFBunchShift);
		locNewEventRFBunch->Set_Content(locInitialEventRFBunch->dTimeSource, locRFTime, locInitialEventRFBunch->dTimeVariance, 0);
		locEventRFBunch = locNewEventRFBunch;
		dRFBunchMap.emplace(locRFBunchShift, locEventRFBunch);
	}
	locParticleCombo->Set_EventRFBunch(locEventRFBunch);
	if(dDebugLevel >= 5)
		cout << "RF Bunch set, shift = " << locRFBunchShift << endl;

	auto locReactionSteps = locReaction->Get_ReactionSteps();
	for(size_t loc_i = 0; loc_i < locReactionSteps.size(); ++loc_i)
	{
		auto locReactionStep = locReactionSteps[loc_i];
		auto locStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(loc_i);
		auto locStepBeamParticle = (loc_i == 0) ? locBeamParticle : nullptr;
		auto locIsProductionVertex = locStepVertexInfo->Get_ProductionVertexFlag();
		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locFullCombo, locStepVertexInfo);
		if(dDebugLevel >= 5)
		{
			cout << "VERTEX PRIMARY COMBO:" << endl;
			Print_SourceCombo(locVertexPrimaryCombo);
		}
		auto locSourceCombo = (loc_i == 0) ? locFullCombo : dSourceComboer->Get_StepSourceCombo(locReaction, loc_i, locVertexPrimaryCombo, locStepVertexInfo->Get_StepIndices().front());
		if(dDebugLevel >= 5)
		{
			cout << "STEP " << loc_i << ", SOURCE COMBO:" << endl;
			Print_SourceCombo(locSourceCombo);
		}

		bool locCreateNeutralErrorMatrixFlag = (locKinFitType != d_NoFit) && ((locKinFitType == d_P4Fit) || !locStepVertexInfo->Get_FittableVertexFlag());
		if(locReactionStep->Get_FinalPIDs(false, d_Neutral, false).empty())
			locCreateNeutralErrorMatrixFlag = false; //no neutrals!

		//reuse step if already created
		auto locStepTuple = std::make_tuple(locSourceCombo, locCreateNeutralErrorMatrixFlag, locIsProductionVertex, locFullCombo, locBeamParticle);
		auto locStepIterator = dComboStepMap.find(locStepTuple);
		if(locStepIterator != dComboStepMap.end())
		{
			locParticleCombo->Add_ParticleComboStep(locStepIterator->second);
			if(dDebugLevel >= 5)
				cout << "step already created, reuse" << endl;
			continue;
		}

		//Create a new step
		auto locParticleComboStep = dResourcePool_ParticleComboStep.Get_Resource();
		locParticleComboStep->Reset();

		//build spacetime vertex
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locVertexPrimaryCombo, locBeamParticle);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locFullCombo, locVertexPrimaryCombo, locBeamParticle);
		auto locPropagatedRFTime = dSourceComboTimeHandler->Calc_PropagatedRFTime(locPrimaryVertexZ, locRFBunchShift, locTimeOffset);
		DLorentzVector locSpacetimeVertex(locVertex, locPropagatedRFTime);
		if(dDebugLevel >= 5)
			cout << "spacetime vertex xyzt: " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

		//Set initial particle
		if((locBeamParticle != nullptr) && (loc_i == 0))
			locParticleComboStep->Set_InitialParticle(locBeamParticle);
		else
			locParticleComboStep->Set_InitialParticle(nullptr);

		//Build final particles
		auto locFinalPIDs = locReactionStep->Get_FinalPIDs();
		map<Particle_t, size_t> locPIDCountMap;
		vector<const DKinematicData*> locFinalParticles;
		for(size_t loc_j = 0; loc_j < locFinalPIDs.size(); ++loc_j)
		{
			//if missing or decaying, nothing to get
			if((int(loc_j) == locReactionStep->Get_MissingParticleIndex()) || (DAnalysis::Get_DecayStepIndex(locReaction, loc_i, loc_j) >= 0))
			{
				locFinalParticles.push_back(nullptr);
				continue;
			}

			//Get source objects, in order
			auto locPID = locFinalPIDs[loc_j];
			auto locPIDIterator = locPIDCountMap.find(locPID);
			if(locPIDIterator != locPIDCountMap.end())
				++(locPIDIterator->second);
			else
				locPIDCountMap.emplace(locPID, 1);

			size_t locPIDCountSoFar = 0;

			auto locSourceParticle = DAnalysis::Get_SourceParticle_ThisStep(locSourceCombo, locPID, locPIDCountMap[locPID], locPIDCountSoFar);
			if(dDebugLevel >= 5)
				cout << "pid index, pid, source particle: " << loc_j << ", " << locPID << ", " << locSourceParticle << endl;

			//build hypo
			if(ParticleCharge(locPID) == 0) //neutral
			{
				auto locNeutralShower = static_cast<const DNeutralShower*>(locSourceParticle);
				auto locHypoTuple = std::make_tuple(locNeutralShower, locPID, locRFBunchShift, locCreateNeutralErrorMatrixFlag, locIsProductionVertex, locFullCombo, locVertexPrimaryCombo, locBeamParticle); //last 4 needed for spacetime vertex
				const DNeutralParticleHypothesis* locNewNeutralHypo = nullptr;

				auto locHypoIterator = dNeutralHypoMap.find(locHypoTuple);
				if(locHypoIterator != dNeutralHypoMap.end())
					locNewNeutralHypo = locHypoIterator->second;
				else
				{
					auto locVertexCovMatrix = locCreateNeutralErrorMatrixFlag ? &dVertexCovMatrix : nullptr;
					locNewNeutralHypo = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralShower, locPID, locEventRFBunch, locSpacetimeVertex, locVertexCovMatrix);
					dNeutralHypoMap.emplace(locHypoTuple, locNewNeutralHypo);
				}

				locFinalParticles.push_back(static_cast<const DKinematicData*>(locNewNeutralHypo));
			}
			else //charged
			{
				auto locChargedTrack = static_cast<const DChargedTrack*>(locSourceParticle);
				auto locHypoTuple = std::make_tuple(locChargedTrack, locPID, locRFBunchShift, locIsProductionVertex, locFullCombo, locVertexPrimaryCombo, locBeamParticle);
				const DChargedTrackHypothesis* locNewChargedHypo = nullptr;

				auto locHypoIterator = dChargedHypoMap.find(locHypoTuple);
				if(locHypoIterator != dChargedHypoMap.end())
					locNewChargedHypo = locHypoIterator->second;
				else
				{
					locNewChargedHypo = Create_ChargedHypo(locChargedTrack, locPID, locPropagatedRFTime, locIsProductionVertex, locVertexPrimaryCombo, locBeamParticle);
					dChargedHypoMap.emplace(locHypoTuple, locNewChargedHypo);
				}

				locFinalParticles.push_back(static_cast<const DKinematicData*>(locNewChargedHypo));
			}
		}
		locParticleComboStep->Set_Contents(locStepBeamParticle, locFinalParticles, locSpacetimeVertex);

		//save it
		locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
		dComboStepMap.emplace(locStepTuple, locParticleComboStep);
	}

	return locParticleCombo;
}

const DChargedTrackHypothesis* DParticleComboCreator::Create_ChargedHypo(const DChargedTrack* locChargedTrack, Particle_t locPID, double locPropagatedRFTime, bool locIsProductionVertex, const DSourceCombo* locVertexPrimaryFullCombo, const DKinematicData* locBeamParticle)
{
	//see if DChargedTrackHypothesis with the desired PID was created by the default factory, AND it passed the PreSelect cuts
	auto locOrigHypo = locChargedTrack->Get_Hypothesis(locPID);
	auto locNewHypo = dChargedTrackHypothesisFactory->Get_Resource();
	locNewHypo->Share_FromInput(locOrigHypo, true, false, true); //share all but timing info

	auto locTrackPOCAX4 = dSourceComboTimeHandler->Get_ChargedParticlePOCAToVertexX4(locOrigHypo, locIsProductionVertex, locVertexPrimaryFullCombo, locBeamParticle);
	locNewHypo->Set_TimeAtPOCAToVertex(locTrackPOCAX4.T());

	locNewHypo->Set_T0(locPropagatedRFTime, locOrigHypo->t0_err(), locOrigHypo->t0_detector());
	locNewHypo->AddAssociatedObject(locChargedTrack);
	dParticleID->Calc_ChargedPIDFOM(locNewHypo);

	return locNewHypo;
}

const DParticleCombo* DParticleComboCreator::Create_KinFitCombo_NewCombo(const DParticleCombo* locOrigCombo, const DReaction* locReaction, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain)
{
	auto locNewCombo = dResourcePool_ParticleCombo.Get_Resource();
	locNewCombo->Reset();
	locNewCombo->Set_KinFitResults(locKinFitResults);
	locNewCombo->Set_EventRFBunch(locOrigCombo->Get_EventRFBunch());
	set<DKinFitParticle*> locOutputKinFitParticles = locKinFitResults->Get_OutputKinFitParticles();

	auto locKinFitType = locKinFitResults->Get_KinFitType();
	for(size_t loc_j = 0; loc_j < locOrigCombo->Get_NumParticleComboSteps(); ++loc_j)
	{
		auto locComboStep = locOrigCombo->Get_ParticleComboStep(loc_j);
		auto locReactionStep = locReaction->Get_ReactionStep(loc_j);

		auto locNewComboStep = dResourcePool_ParticleComboStep.Get_Resource();
		locNewComboStep->Reset();
		locNewCombo->Add_ParticleComboStep(locNewComboStep);

		locNewComboStep->Set_MeasuredParticleComboStep(locComboStep);
		locNewComboStep->Set_SpacetimeVertex(locComboStep->Get_SpacetimeVertex()); //overridden if kinematic fit

		//INITIAL PARTICLE
		auto locInitialParticle_Measured = locComboStep->Get_InitialParticle_Measured();
		if(locInitialParticle_Measured != nullptr) //set beam photon
		{
			auto locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locInitialParticle_Measured);
			if(locKinFitParticle == NULL) //not used in kinfit!!
				locNewComboStep->Set_InitialParticle(locInitialParticle_Measured);
			else //create a new one
				locNewComboStep->Set_InitialParticle(Create_BeamPhoton_KinFit(static_cast<const DBeamPhoton*>(locInitialParticle_Measured), locKinFitParticle));
		}
		else //decaying particle! //set here for initial state, and in previous step for final state
			Set_DecayingParticles(locReaction, locNewCombo, locOrigCombo, loc_j, locNewComboStep, locKinFitChain, locKinFitResults);

		//FINAL PARTICLES
		for(size_t loc_k = 0; loc_k < locComboStep->Get_NumFinalParticles(); ++loc_k)
		{
			auto locKinematicData_Measured = locComboStep->Get_FinalParticle_Measured(loc_k);
			if(locReactionStep->Get_MissingParticleIndex() == int(loc_k)) //missing!
			{
				set<DKinFitParticle*> locMissingParticles = locKinFitResults->Get_OutputKinFitParticles(d_MissingParticle);
				if(!locMissingParticles.empty())
				{
					DKinematicData* locNewKinematicData = Build_KinematicData(*locMissingParticles.begin(), locKinFitType, locComboStep->Get_SpacetimeVertex().Vect());
					locNewComboStep->Add_FinalParticle(locNewKinematicData);
				}
				else //not used in kinfit: do not create: NULL
					locNewComboStep->Add_FinalParticle(NULL);
			}
			else if(locKinematicData_Measured == nullptr) //decaying
				locNewComboStep->Add_FinalParticle(NULL); //is set later, when it's in the initial state
			else if(ParticleCharge(locKinematicData_Measured->PID()) == 0) //neutral
			{
				auto locNeutralHypo = static_cast<const DNeutralParticleHypothesis*>(locKinematicData_Measured);
				//might have used neutral shower OR neutral hypo. try hypo first
				auto locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locNeutralHypo);
				if(locKinFitParticle == NULL)
					locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locNeutralHypo->Get_NeutralShower());
				if(locKinFitParticle == NULL) //not used in kinfit!!
					locNewComboStep->Add_FinalParticle(locKinematicData_Measured);
				else //create a new one
					locNewComboStep->Add_FinalParticle(Create_NeutralHypo_KinFit(locNeutralHypo, locKinFitParticle));
			}
			else //charged
			{
				auto locChargedHypo = static_cast<const DChargedTrackHypothesis*>(locKinematicData_Measured);
				auto locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locChargedHypo);
				if(locKinFitParticle == NULL) //not used in kinfit!!
					locNewComboStep->Add_FinalParticle(locKinematicData_Measured);
				else //create a new one
				{
					auto locChargedTrack = static_cast<const DChargedTrack*>(locComboStep->Get_FinalParticle_SourceObject(loc_k));
					auto locNewHypo = Create_ChargedHypo_KinFit(locChargedTrack, locKinematicData_Measured->PID(), locKinFitParticle);
					locNewComboStep->Add_FinalParticle(locNewHypo);
				}
			}
		}

		//SPACETIME VERTEX
		Set_SpacetimeVertex(locReaction, locNewCombo, locNewComboStep, loc_j, locKinFitResults, locKinFitChain);
	}

	return locNewCombo;
}

void DParticleComboCreator::Set_DecayingParticles(const DReaction* locReaction, const DParticleCombo* locNewParticleCombo, const DParticleCombo* locOldParticleCombo, size_t locStepIndex, DParticleComboStep* locNewParticleComboStep, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults)
{
	DKinFitParticle* locKinFitParticle = Get_DecayingParticle(locReaction, locOldParticleCombo, locStepIndex, locKinFitChain, locKinFitResults);
	if(locKinFitParticle == NULL) //not used in fit
	{
		locNewParticleComboStep->Set_InitialParticle(NULL);
		return; //no need to back-set NULL: was set to NULL by default
	}

	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	DVector3 locEventVertex = locOldParticleCombo->Get_ParticleComboStep(0)->Get_SpacetimeVertex().Vect();
	Particle_t locPID = PDGtoPType(locKinFitParticle->Get_PID());
	auto locFromStepIndex = DAnalysis::Get_InitialParticleDecayFromIndices(locReaction, locStepIndex).first;

	DKinematicData* locKinematicData_Common = Build_KinematicData(locKinFitParticle, locKinFitType, locEventVertex);
	bool locCreate2ndObjectFlag = (IsDetachedVertex(locPID) && (locStepIndex != 0) && (locFromStepIndex >= 0));
	DKinematicData* locKinematicData_Position = locCreate2ndObjectFlag ? Build_KinematicData(locKinFitParticle, locKinFitType, locEventVertex) : locKinematicData_Common;
	if(locKinFitParticle->Get_CommonVxParamIndex() >= 0)
		dKinFitUtils->Propagate_TrackInfoToCommonVertex(locKinematicData_Common, locKinFitParticle, &locKinFitResults->Get_VXi());

	bool locAtProdVertexFlag = locKinFitParticle->Get_VertexP4AtProductionVertex();
	DKinematicData* locKinematicData_InitState = locAtProdVertexFlag ? locKinematicData_Common : locKinematicData_Position;
	DKinematicData* locKinematicData_FinalState = locAtProdVertexFlag ? locKinematicData_Position : locKinematicData_Common;

	locNewParticleComboStep->Set_InitialParticle(locKinematicData_InitState);

	//now, back-set the particle at the other vertex
	if((locStepIndex == 0) || (locFromStepIndex < 0))
		return; //no other place to set it

	DParticleComboStep* locParticleComboStep = const_cast<DParticleComboStep*>(locNewParticleCombo->Get_ParticleComboStep(locFromStepIndex));
	//find where it is the decaying particle
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		auto locDecayStepIndex = DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i);
		if(locDecayStepIndex != int(locStepIndex))
			continue;
		locParticleComboStep->Set_FinalParticle(locKinematicData_FinalState, loc_i);
		break;
	}
}

DKinFitParticle* DParticleComboCreator::Get_DecayingParticle(const DReaction* locReaction, const DParticleCombo* locOldParticleCombo, size_t locComboStepIndex, const DKinFitChain* locKinFitChain, const DKinFitResults* locKinFitResults)
{
	auto locReactionStep = locReaction->Get_ReactionStep(locComboStepIndex);
	Particle_t locPID = locReactionStep->Get_InitialPID();
	if(!IsFixedMass(locPID))
		return NULL;

	//find which step in the DKinFitChain this combo step corresponds to
	for(size_t loc_i = 0; loc_i < locKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
	{
		const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(loc_i);

		//loop over init particles to get the decaying particle (if present)
		DKinFitParticle* locDecayingParticle = NULL;
		auto locInitialParticles = locKinFitChainStep->Get_InitialParticles();
		for(auto locKinFitParticle : locInitialParticles)
		{
			if(locKinFitParticle == nullptr)
				continue;
			if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
				continue; //not a decaying particle
			locDecayingParticle = locKinFitParticle;
			break;
		}
		if(locDecayingParticle == NULL)
			continue; //no decaying particles in this step

		if(PDGtoPType(locDecayingParticle->Get_PID()) != locPID)
			continue; //wrong PID

		//may still not be correct particle. compare decay products: if any of the combo step particles are in the kinfit step, this is it
			//if all step final particles are decaying, then dive down: through steps
			//if any step decay product at any step is located as any decay product at any step in the kinfit chain: then matches

		//get all measured products, then just pick the first one to search for
		auto locMeasuredParticles = locOldParticleCombo->Get_DecayChainParticles_Measured(locReaction, locComboStepIndex);
		const DKinematicData* locMeasuredParticle = locMeasuredParticles[0];
		DKinFitParticle* locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locMeasuredParticle);
		if(locKinFitParticle == NULL) //null: neutral shower. Use shower object
		{
			const DNeutralShower* locNeutralShower = static_cast<const DNeutralParticleHypothesis*>(locMeasuredParticle)->Get_NeutralShower();
			locKinFitParticle = locKinFitResults->Get_OutputKinFitParticle(locNeutralShower);
		}

		if(!Search_ForParticleInDecay(locKinFitChain, loc_i, locKinFitParticle))
			continue;

		//Found!
		return locDecayingParticle;
	}

	return NULL;
}

bool DParticleComboCreator::Search_ForParticleInDecay(const DKinFitChain* locKinFitChain, size_t locStepToSearch, DKinFitParticle* locParticleToFind)
{
	const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(locStepToSearch);
	auto locFinalParticles = locKinFitChainStep->Get_FinalParticles();
	std::sort(locFinalParticles.begin(), locFinalParticles.end());
	if(std::binary_search(locFinalParticles.begin(), locFinalParticles.end(), locParticleToFind))
		return true; //found it

	//else loop over final state, diving through decays
	for(auto locKinFitParticle : locFinalParticles)
	{
		if(locKinFitParticle == nullptr)
			continue;
		if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //not a decaying particle

		int locDecayStepIndex = locKinFitChain->Get_DecayStepIndex(locKinFitParticle);
		if(Search_ForParticleInDecay(locKinFitChain, locDecayStepIndex, locParticleToFind))
			return true; //found in in subsequent step
	}
	return false; //not found (yet)
}

void DParticleComboCreator::Set_SpacetimeVertex(const DReaction* locReaction, const DParticleCombo* locNewParticleCombo, DParticleComboStep* locNewParticleComboStep, size_t locStepIndex, const DKinFitResults* locKinFitResults, const DKinFitChain* locKinFitChain) const
{
	DKinFitType locKinFitType = locReaction->Get_KinFitType();
	if((locKinFitType == d_NoFit) || (locKinFitType == d_P4Fit))
		return; //neither vertex nor time was fit: no update to give

	//Position & Time
	const DKinematicData* locKinematicData = locNewParticleComboStep->Get_InitialParticle();
	if(locKinematicData != NULL)
	{
		DLorentzVector locSpacetimeVertex(locKinematicData->position(), locKinematicData->time());
		locNewParticleComboStep->Set_SpacetimeVertex(locSpacetimeVertex);
		return;
	}

	//mass not fixed: if not initial step, get spacetime vertex from step where this particle was produced
	if(locStepIndex != 0)
	{
		//get spacetime vertex from step where this particle was produced
		auto locDecayFromStepIndex = Get_InitialParticleDecayFromIndices(locReaction, locStepIndex).first;
		const DParticleComboStep* locPreviousParticleComboStep = locNewParticleCombo->Get_ParticleComboStep(locDecayFromStepIndex);
		locNewParticleComboStep->Set_SpacetimeVertex(locPreviousParticleComboStep->Get_SpacetimeVertex());
		return;
	}

	//instead, get from common vertex of the other particles
	const DKinFitParticle* locFinalKinFitParticle = nullptr;
	auto locAllParticles = locKinFitChain->Get_AllParticles();
	for(auto& locFinalKinFitParticle : locAllParticles)
	{
		if(locFinalKinFitParticle != nullptr)
			break;
	}
	if(locFinalKinFitParticle == nullptr)
		return;

	//need the spacetime vertex at the production vertex of the particle grabbed
	TLorentzVector locSpacetimeVertex;
	if(locFinalKinFitParticle->Get_VertexP4AtProductionVertex()) //"position" is at production vertex
		locSpacetimeVertex = locFinalKinFitParticle->Get_SpacetimeVertex();
	else //"position" is at decay vertex
		locSpacetimeVertex = locFinalKinFitParticle->Get_CommonSpacetimeVertex(); //get production
	locNewParticleComboStep->Set_SpacetimeVertex(locSpacetimeVertex);
}

const DBeamPhoton* DParticleComboCreator::Create_BeamPhoton_KinFit(const DBeamPhoton* locBeamPhoton, const DKinFitParticle* locKinFitParticle)
{
	auto locBeamIterator = dKinFitBeamPhotonMap.find(locKinFitParticle);
	if(locBeamIterator != dKinFitBeamPhotonMap.end())
		return locBeamIterator->second;

	DBeamPhoton* locNewBeamPhoton = dBeamPhotonfactory->Get_Resource();
	dKinFitBeamPhotonMap.emplace(locKinFitParticle, locNewBeamPhoton);

	locNewBeamPhoton->dCounter = locBeamPhoton->dCounter;
	locNewBeamPhoton->dSystem = locBeamPhoton->dSystem;
	locNewBeamPhoton->setPID(Gamma);
	locNewBeamPhoton->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewBeamPhoton->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewBeamPhoton->setTime(locKinFitParticle->Get_Time());
	locNewBeamPhoton->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());
	return locNewBeamPhoton;
}

const DChargedTrackHypothesis* DParticleComboCreator::Create_ChargedHypo_KinFit(const DChargedTrack* locChargedTrack, Particle_t locPID, const DKinFitParticle* locKinFitParticle)
{
	auto locOrigHypo = locChargedTrack->Get_Hypothesis(locPID);
	auto locHypoIterator = dKinFitChargedHypoMap.find(locKinFitParticle);
	if(locHypoIterator != dKinFitChargedHypoMap.end())
		return locHypoIterator->second;

	//even if vertex is not fit, p4 is fit: different beta: update time info
	auto locNewHypo = dChargedTrackHypothesisFactory->Get_Resource();
	dKinFitChargedHypoMap.emplace(locKinFitParticle, locNewHypo);
	locNewHypo->setPID(locPID);

	//p3 & v3
	TVector3 locFitMomentum = locKinFitParticle->Get_Momentum();
	TVector3 locFitVertex = locKinFitParticle->Get_Position();
	locNewHypo->setMomentum(DVector3(locFitMomentum.X(), locFitMomentum.Y(), locFitMomentum.Z()));
	locNewHypo->setPosition(DVector3(locFitVertex.X(), locFitVertex.Y(), locFitVertex.Z()));

	//t & error matrix
	locNewHypo->setTime(locKinFitParticle->Get_Time());
	locNewHypo->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());

	//if timing was kinfit, there is no chisq to calculate: forced to correct time
	//therefore, just use measured timing info (pre-kinfit)
	if(locKinFitParticle->Get_CommonTParamIndex() >= 0)
	{
		locNewHypo->Share_FromInput(locOrigHypo, true, true, false); //share all but kinematics
		return locNewHypo;
	}

	//only share tracking info (not timing or kinematics)
	locNewHypo->Share_FromInput(locOrigHypo, true, false, false);

	//update timing info
	if(locKinFitParticle->Get_CommonVxParamIndex() >= 0) //a vertex was fit
	{
		//all kinematics propagated to vertex position
		auto locPropagatedRFTime = locOrigHypo->t0() + (locFitVertex.Z() - locOrigHypo->position().Z())/SPEED_OF_LIGHT;
		locNewHypo->Set_T0(locPropagatedRFTime, locOrigHypo->t0_err(), locOrigHypo->t0_detector());
		locNewHypo->Set_TimeAtPOCAToVertex(locNewHypo->time());
	}
	else //only momentum has changed (and thus time)
	{
		locNewHypo->Set_T0(locOrigHypo->t0(), locOrigHypo->t0_err(), locOrigHypo->t0_detector());
		locNewHypo->Set_TimeAtPOCAToVertex(locOrigHypo->Get_TimeAtPOCAToVertex() + locNewHypo->time() - locOrigHypo->time());
	}

	locNewHypo->AddAssociatedObject(locChargedTrack);
	dParticleID->Calc_ChargedPIDFOM(locNewHypo);
	return locNewHypo;
}

const DNeutralParticleHypothesis* DParticleComboCreator::Create_NeutralHypo_KinFit(const DNeutralParticleHypothesis* locOrigHypo, DKinFitParticle* locKinFitParticle)
{
	auto locHypoIterator = dKinFitNeutralHypoMap.find(locKinFitParticle);
	if(locHypoIterator != dKinFitNeutralHypoMap.end())
		return locHypoIterator->second;

	auto locNewHypo = dNeutralParticleHypothesisFactory->Get_Resource();
	dKinFitNeutralHypoMap.emplace(locKinFitParticle, locNewHypo);
	locNewHypo->setPID(locOrigHypo->PID());
	locNewHypo->Set_NeutralShower(locOrigHypo->Get_NeutralShower());

	//p3 & v3
	TVector3 locFitMomentum = locKinFitParticle->Get_Momentum();
	TVector3 locFitVertex = locKinFitParticle->Get_Position();
	locNewHypo->setMomentum(DVector3(locFitMomentum.X(), locFitMomentum.Y(), locFitMomentum.Z()));
	locNewHypo->setPosition(DVector3(locFitVertex.X(), locFitVertex.Y(), locFitVertex.Z()));

	//t & error matrix
	locNewHypo->setTime(locKinFitParticle->Get_Time());
	locNewHypo->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());

	//if timing was kinfit, there is no chisq to calculate: forced to correct time
	//also, if vertex was not kinfit, no chisq to calc either: photon beta is still 1, no chisq for massive neutrals
	//therefore, just use measured timing info (pre-kinfit)
	if((locKinFitParticle->Get_CommonVxParamIndex() < 0) || (locKinFitParticle->Get_CommonTParamIndex() >= 0))
	{
		locNewHypo->Share_FromInput(locOrigHypo, true, false); //share timing but not kinematics
		return locNewHypo;
	}

	//update timing info
	auto locPropagatedRFTime = locOrigHypo->t0() + (locFitVertex.Z() - locOrigHypo->position().Z())/SPEED_OF_LIGHT;
	locNewHypo->Set_T0(locPropagatedRFTime, locOrigHypo->t0_err(), locOrigHypo->t0_detector());

	// Calculate DNeutralParticleHypothesis FOM
	unsigned int locNDF = 0;
	double locChiSq = 0.0;
	double locFOM = -1.0; //undefined for non-photons
	if(locNewHypo->PID() == Gamma)
	{
		double locTimePull = 0.0;
		//for this calc: if rf time part of timing constraint, don't use locKinFitParticle->Get_CommonTime() for chisq calc!!!
		locChiSq = dParticleID->Calc_TimingChiSq(locNewHypo, locNDF, locTimePull);
		locFOM = TMath::Prob(locChiSq, locNDF);
	}
	locNewHypo->Set_ChiSq_Overall(locChiSq, locNDF, locFOM);

	return locNewHypo;
}

DKinematicData* DParticleComboCreator::Build_KinematicData(DKinFitParticle* locKinFitParticle, DKinFitType locKinFitType, DVector3 locPreKinFitVertex)
{
	DKinematicData* locKinematicData = dResourcePool_KinematicData.Get_Resource();
	locKinematicData->Reset();
	locKinematicData->setPID(PDGtoPType(locKinFitParticle->Get_PID()));
	locKinematicData->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	if((locKinFitType == d_P4Fit) || (locKinFitType == d_NoFit))
		locKinematicData->setPosition(locPreKinFitVertex);
	else
		locKinematicData->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locKinematicData->setTime(locKinFitParticle->Get_Time());
	if(locKinFitParticle->Get_CovarianceMatrix() != NULL)
		locKinematicData->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());

	return locKinematicData;
}

const DParticleCombo* DParticleComboCreator::Build_ThrownCombo(JEventLoop* locEventLoop)
{
	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	if(dAnalysisUtilities == nullptr)
		locEventLoop->GetSingle(dAnalysisUtilities);
	dAnalysisUtilities->Get_ThrownParticleSteps(locEventLoop, locThrownSteps);
	if(locThrownSteps.empty())
		return nullptr;

 	vector<const DReaction*> locReactions;
	locEventLoop->Get(locReactions, "Thrown");
	return Build_ThrownCombo(locEventLoop, locReactions[0], locThrownSteps);
}

const DParticleCombo* DParticleComboCreator::Build_ThrownCombo(JEventLoop* locEventLoop, const DReaction* locThrownReaction, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps)
{
	auto locThrownComboTuple = std::make_tuple((const DReactionVertexInfo*)nullptr, (const DSourceCombo*)nullptr, (const DKinematicData*)nullptr, 1, true);
	auto locComboMapIterator = dComboMap.find(locThrownComboTuple);
	if(locComboMapIterator != dComboMap.end())
		return locComboMapIterator->second;

 	vector<const DMCReaction*> locMCReactions;
	locEventLoop->Get(locMCReactions);

 	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Thrown");

	DParticleCombo* locParticleCombo = dResourcePool_ParticleCombo.Get_Resource();
	DParticleComboStep* locParticleComboStep = dResourcePool_ParticleComboStep.Get_Resource();
	locParticleCombo->Set_EventRFBunch(locEventRFBunches[0]);

	locParticleComboStep->Set_InitialParticle(&locMCReactions[0]->beam);
	locParticleComboStep->Set_SpacetimeVertex(DLorentzVector(locMCReactions[0]->beam.position(), locMCReactions[0]->beam.time()));

	for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
	{
		if(loc_i != 0) //else beam & target already set
		{
			auto locMCThrown = locThrownSteps[loc_i].first;
			locParticleComboStep = dResourcePool_ParticleComboStep.Get_Resource();

			int locInitialParticleDecayFromStepIndex = -1; //the step where this particle is produced at
			for(size_t loc_j = 0; loc_j < loc_i; ++loc_j)
			{
				for(size_t loc_k = 0; loc_k < locThrownSteps[loc_j].second.size(); ++loc_k)
				{
					if(locMCThrown != locThrownSteps[loc_j].second[loc_k])
						continue;
					locInitialParticleDecayFromStepIndex = loc_j;
					break;
				}
				if(locInitialParticleDecayFromStepIndex != -1)
					break;
			}
			locParticleComboStep->Set_InitialParticle(locMCThrown);
			locParticleComboStep->Set_SpacetimeVertex(DLorentzVector(locMCThrown->position(), locMCThrown->time()));
		}
		for(size_t loc_j = 0; loc_j < locThrownSteps[loc_i].second.size(); ++loc_j)
			locParticleComboStep->Add_FinalParticle(locThrownSteps[loc_i].second[loc_j]);
		locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
	}

	dComboMap.emplace(locThrownComboTuple, locParticleCombo);
	return locParticleCombo;
}

}
