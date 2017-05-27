#include "ANALYSIS/DParticleComboCreator.h"
#include "ANALYSIS/DSourceComboer.h"

namespace DAnalysis
{

void DParticleComboCreator::DParticleComboCreator(JEventLoop* locEventLoop, const DSourceComboer* locSourceComboer, const DSourceComboTimeHandler* locSourceComboTimeHandler, const DSourceComboVertexer* dSourceComboVertexer) :
		dSourceComboer(locSourceComboer), dSourceComboTimeHandler(locSourceComboTimeHandler), dSourceComboVertexer(locSourceComboVertexer)
{
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
	dResourcePool_EventRFBunch.clear();

	for(const auto& locStepPair : dComboStepMap)
		dResourcePool_ParticleComboStep.Recycle(locStepPair.second);
	dComboStepMap.clear();

	for(const auto& locComboPair : dComboMap)
		dResourcePool_ParticleCombo.Recycle(locComboPair.second);
	dComboMap.clear();

	for(const auto& locHypoPair : dChargedHypoMap)
		dChargedTrackHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dChargedHypoMap.clear();

	for(const auto& locHypoPair : dNeutralHypoMap)
		dNeutralParticleHypothesisFactory->Recycle_Hypothesis(locHypoPair.second);
	dNeutralHypoMap.clear();
}

bool DParticleComboCreator::Get_CreateNeutralErrorMatrixFlag_Combo(const DReactionVertexInfo* locReactionVertexInfo, DKinFitType locKinFitType)
{
	//is there at least one dangling vertex that has neutrals?
	auto locDanglingIterator = dDanglingNeutralsFlagMap.find(locReactionVertexInfo);
	bool locDanglingNeutralsFlag = false;
	if(locDanglingIterator != dDanglingNeutralsFlagMap.end())
		locDanglingNeutralsFlag = locDanglingIterator->second;
	else
	{
		for(auto locStepVertexInfo : locReactionVertexInfo->Get_StepVertexInfos())
		{
			if(!locStepVertexInfo->Get_DanglingVertexFlag())
				continue;
			if(!locStepVertexInfo->Get_OnlyConstrainTimeParticles().empty())
			{
				locDanglingNeutralsFlag = true; //photons
				break;
			}
			if(!locStepVertexInfo->Get_NoConstrainParticles(d_FinalState, d_Neutral, false, false, false).empty())
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
	auto locCreateNeutralErrorMatrixFlag_Combo = Get_CreateNeutralErrorMatrixFlag_Combo(locReactionVertexInfo, locKinFitType);
	auto locComboTuple = std::make_tuple(locReactionVertexInfo, locFullCombo, locBeamParticle, locRFBunchShift, locCreateNeutralErrorMatrixFlag_Combo);
	auto locComboIterator = dComboMap.find(locComboTuple);
	if(locComboIterator != dComboMap.end())
		return locComboIterator->second;

	auto locParticleCombo = dResourcePool_ParticleCombo.Get_Resource();
	locParticleCombo->Reset();
	dComboMap.emplace(locComboTuple, locParticleCombo);

	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();
	auto locPrimaryVertexZ = dSourceComboVertexer->Get_PrimaryVertex(locReactionVertexInfo, locFullCombo, locBeamParticle).Z();

	//Get/Create RF Bunch
	const DEventRFBunch* locEventRFBunch = nullptr;
	auto locRFIterator = dRFBunchMap.find(locRFBunchShift);
	if(locRFIterator != dRFBunchMap.end())
		locEventRFBunch = locRFIterator->second;
	else
	{
		auto locInitialEventRFBunch = dSourceComboTimeHandler->Get_InitialEventRFBunch();
		auto locEventRFBunch = dResourcePool_EventRFBunch.Get_Resource();
		auto locRFTime = dSourceComboTimeHandler->Calc_RFTime(locRFBunchShift);
		locEventRFBunch->Set_Content(locInitialEventRFBunch->dTimeSource, locRFTime, locInitialEventRFBunch->dTimeVariance, 0);
	}
	locParticleCombo->Set_EventRFBunch(locEventRFBunch);

	auto locReactionSteps = locReaction->Get_ReactionSteps();
	auto locPreviousStepSourceCombo = locFullCombo;
	for(size_t loc_i = 0; loc_i < locReactionSteps.size(); ++loc_i)
	{
		auto locReactionStep = locReactionSteps[loc_i];
		auto locStepBeamParticle = (loc_i == 0) ? locBeamParticle : nullptr;
		auto locIsProductionVertex = (loc_i == 0) ? locIsPrimaryProductionVertex : false;
		auto locSourceCombo = (loc_i == 0) ? locFullCombo : dSourceComboer->Get_StepSourceCombo(locReaction, loc_i, locPreviousStepSourceCombo, loc_i - 1);

		auto locStepVertexInfo = locReactionVertexInfo->Get_StepVertexInfo(loc_i);
		auto locVertexPrimaryCombo = dSourceComboer->Get_VertexPrimaryCombo(locFullCombo, locStepVertexInfo);
		bool locCreateNeutralErrorMatrixFlag = (locKinFitType != d_NoFit) && ((locKinFitType == d_P4Fit) || locStepVertexInfo->Get_DanglingVertexFlag());
		if(locReactionStep->Get_FinalPIDs(false, d_Neutral, false).empty())
			locCreateNeutralErrorMatrixFlag = false; //no neutrals!

		//reuse step if already created
		auto locStepTuple = std::make_tuple(locSourceCombo, locCreateNeutralErrorMatrixFlag, locIsProductionVertex, locFullCombo, locBeamParticle);
		auto locStepIterator = dComboStepMap.find(locStepTuple);
		if(locStepIterator != dComboStepMap.end())
		{
			locParticleCombo->Add_ParticleComboStep(locStepIterator->second);
			locPreviousStepSourceCombo = locSourceCombo;
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

		//Build final particles
		auto locFinalPIDs = locReactionStep->Get_FinalPIDs();
		unordered_map<Particle_t, size_t> locPIDCountMap;
		vector<const DKinematicData*> locFinalParticles;
		for(size_t loc_j = 0; loc_j < locFinalPIDs.size(); ++loc_j)
		{
			if(loc_j == locReactionStep->Get_MissingParticleIndex())
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

			//build hypo
			if(ParticleCharge(locPID) == 0)
			{
				auto locNeutralShower = static_cast<const DNeutralShower*>(locSourceParticle);
				auto locHypoTuple = std::make_tuple(locNeutralShower, locPID, locRFBunchShift, locCreateNeutralErrorMatrixFlag, locIsProductionVertex, locFullCombo, locVertexPrimaryCombo, locBeamParticle); //last 4 needed for spacetime vertex
				const DNeutralParticleHypothesis* locNewNeutralHypo = nullptr;

				auto locHypoIterator = dNeutralHypoMap.find(locHypoTuple);
				if(locHypoIterator != dNeutralHypoMap.end())
					locNewNeutralHypo = locHypoIterator->second;
				else
				{
					auto locVertexCovMatrix = locCreateNeutralErrorMatrixFlag ? *dVertexCovMatrix : nullptr;
					locNewNeutralHypo = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralShower, locPID, locEventRFBunch, locSpacetimeVertex, locVertexCovMatrix);
					dNeutralHypoMap.emplace(locHypoTuple, locNewNeutralHypo);
				}

				locFinalParticles.push_back(static_cast<DKinematicData*>(locNewNeutralHypo));
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

				locFinalParticles.push_back(static_cast<DKinematicData*>(locNewChargedHypo));
			}
		}

		locParticleComboStep->Set_Contents(locStepBeamParticle, locFinalParticles, locSpacetimeVertex);
		if(loc_i == 0)
			locParticleComboStep->Set_InitialParticle(locBeamParticle);

		//save it
		locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
		dComboStepMap.emplace(locStepTuple, locParticleComboStep);

		locPreviousStepSourceCombo = locSourceCombo;
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
	dParticleID->Calc_ChargedPIDFOM(locNewHypo);

	return locNewHypo;
}

const DBeamPhoton* DParticleComboCreator::Build_BeamPhoton_KinFit(const DBeamPhoton* locBeamPhoton, DKinFitParticle* locKinFitParticle)
{
	DBeamPhoton* locNewBeamPhoton = dBeamPhotonfactory->Get_Resource();
	locNewBeamPhoton->dCounter = locBeamPhoton->dCounter;
	locNewBeamPhoton->dSystem = locBeamPhoton->dSystem;
	locNewBeamPhoton->setMomentum(DVector3(locKinFitParticle->Get_Momentum().X(),locKinFitParticle->Get_Momentum().Y(),locKinFitParticle->Get_Momentum().Z()));
	locNewBeamPhoton->setPosition(DVector3(locKinFitParticle->Get_Position().X(),locKinFitParticle->Get_Position().Y(),locKinFitParticle->Get_Position().Z()));
	locNewBeamPhoton->setTime(locKinFitParticle->Get_Time());
	locNewBeamPhoton->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());
	return locNewBeamPhoton;
}



DChargedTrackHypothesis* DParticleComboCreator::Build_ChargedTrackHypothesis(const DChargedTrackHypothesis* locOrigHypo, DKinFitParticle* locKinFitParticle, bool locVertexFitFlag)
{
	auto locNewHypo = dChargedTrackHypothesisFactory->Get_Resource();

	//p3 & v3
	TVector3 locFitMomentum = locKinFitParticle->Get_Momentum();
	TVector3 locFitVertex = locKinFitParticle->Get_Position();
	locNewHypo->setMomentum(DVector3(locFitMomentum.X(), locFitMomentum.Y(), locFitMomentum.Z()));
	locNewHypo->setPosition(DVector3(locFitVertex.X(), locFitVertex.Y(), locFitVertex.Z()));

	//t & error matrix
	locNewHypo->setTime(locKinFitParticle->Get_Time());
	locNewHypo->setErrorMatrix(locKinFitParticle->Get_CovarianceMatrix());

	if(!locVertexFitFlag)
		locNewHypo->Share_FromInput(locOrigHypo, true, true, false); //share all but kinematics
	else //new timing info
	{
//OK
		locNewHypo->Share_FromInput(locOrigHypo, true, false, false); //only share tracking info (not timing or kinematics)
		auto locPropagatedRFTime = locOrigHypo->t0() + (locCommonVertex.Z() - locOrigHypo->position().Z())/SPEED_OF_LIGHT;
		locNewHypo->Set_T0(locPropagatedRFTime, locOrigHypo->t0_err(), locOrigHypo->t0_detector());

		//BE CAREFUL HERE!
		if(locSpacetimeFitFlag)
		{

		}
		else
		{
			auto locDeltaPathLength = ;
			auto locTimeAtPOCAToVertex = locKinFitParticle->Get_Time();
			locNewHypo->Set_TimeAtPOCAToVertex(locTimeAtPOCAToVertex);

		}

		//for this calc: if rf time part of timing constraint, don't use locKinFitParticle->Get_Time() for chisq calc!!!
		dParticleID->Calc_ChargedPIDFOM(locNewHypo);
	}

	return locNewHypo;
}

}

#endif // DParticleComboCreator_h
