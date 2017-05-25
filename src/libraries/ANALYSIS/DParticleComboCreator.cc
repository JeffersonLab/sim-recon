#include "ANALYSIS/DParticleComboCreator.h"

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
}


const DParticleCombo* DParticleComboCreator::Build_ParticleCombo(JEventLoop* locEventLoop, const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift)
{
	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locPreviousStepSourceCombo = locFullCombo;
	auto locParticleCombo = dResourcePool_ParticleCombo.Get_Resource();
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

	auto locReactionSteps = locReaction->Get_ReactionSteps();
	for(size_t loc_i = 0; loc_i < locReactionSteps.size(); ++loc_i)
	{
		auto locReactionStep = locReactionSteps[loc_i];
		auto locStepBeamParticle = (loc_i == 0) ? locBeamParticle : nullptr;
		auto locIsProductionVertex = (loc_i == 0) ? locIsPrimaryProductionVertex : false;
		auto locSourceCombo = (loc_i == 0) ? locFullCombo : dSourceComboer->Get_StepSourceCombo(locReaction, loc_i, locPreviousStepSourceCombo, loc_i - 1);
		auto locStepTuple = std::make_tuple(locReactionStep, locSourceCombo, locStepBeamParticle);

		//reuse step if already created
		auto locStepIterator = dComboStepMap.find(locStepTuple);
		if(locStepIterator != dComboStepMap.end())
		{
			locParticleCombo->Add_ParticleComboStep(locStepIterator->second);
			continue;
		}

		//Create a new step
		auto locParticleComboStep = dResourcePool_ParticleComboStep.Get_Resource();
		locParticleComboStep->Reset();

		//build spacetime vertex
		auto locVertex = dSourceComboVertexer->Get_Vertex(locIsProductionVertex, locSourceCombo, locBeamParticle);
		auto locTimeOffset = dSourceComboVertexer->Get_TimeOffset(locIsProductionVertex, locFullCombo, locSourceCombo, locBeamParticle);
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
				auto locNeutralShower = static_cast<const DNeutralShower*>(locNeutralShower);
				auto locHypoTuple = std::make_tuple(locNeutralShower, locPID, locRFBunchShift, locIsProductionVertex, locFullCombo, locSourceCombo, locBeamParticle); //last 4 needed for spacetime vertex
	//use tuple, recycle resources when done!
				auto locNeutralHypo = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locEventLoop, locNeutralShower, locPID, locEventRFBunch, locSpacetimeVertex, locVertexCovMatrix);
				locFinalParticles.push_back(static_cast<DKinematicData*>(locNeutralHypo));
			}
			else //charged
			{

			}
		}

		locParticleComboStep->Set_Contents(locReactionStep, locStepBeamParticle, locFinalParticles, locSpacetimeVertex);
		if(loc_i == 0)
			locParticleComboStep->Set_InitialParticle(locBeamParticle);

		// SET FINAL PARTICLES:
		inline void Add_FinalParticle();

		//save it
		locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
		dComboStepMap.emplace(locStepTuple, locParticleComboStep);

		locPreviousStepSourceCombo = locSourceCombo;
	}
}

const DChargedTrackHypothesis* DParticleComboCreator::Create_Charged(const DChargedTrack* locChargedTrack, Particle_t locPID, const DEventRFBunch* locEventRFBunch, double locPropagatedRFTime)
{
	//see if DChargedTrackHypothesis with the desired PID was created by the default factory, AND it passed the PreSelect cuts
	auto locOrigHypo = locChargedTrack->Get_Hypothesis(locPID);
	auto locNewHypo = dChargedTrackHypothesisFactory->Get_Resource();
	locNewHypo->Share_FromInput(locOrigHypo, true, false, true); //share all but timing info
	double locTrackPOCATime = locEventRFBunch->
	locNewHypo->Set_T0(locPropagatedRFTime, locOrigHypo->t0_err(), locOrigHypo->t0_detector());

//dTimeAtPOCAToVertex: what about DSelector? Need to save position as well? eh?
	//So this means that the delta-t used for the initial timing cuts is NOT available later!!

	unsigned int locTimingNDF = 0;
	double locTimingPull = 0.0;
	double locTimingChiSq = 0.0;
	if((locNewHypo->t0_detector() != SYS_NULL) && (locNewHypo->t1_detector() != SYS_NULL))
	{
		double locDeltaT = ;
		double locStartTimeError = locChargedHypo->t0_err();
		double locTimeDifferenceVariance = (*locChargedHypo->errorMatrix())(6, 6) + locStartTimeError*locStartTimeError;
		locTimingPull = (locChargedHypo->t0() - locChargedHypo->time())/sqrt(locTimeDifferenceVariance);
		locTimingNDF = 1;
		locTimingChiSq = locTimingPull*locTimingPull;
	}


	unsigned int locNDF_Total = locChargedTrackHypothesis->Get_NDF_Timing() + locChargedTrackHypothesis->Get_NDF_DCdEdx();
	double locChiSq_Total = locChargedTrackHypothesis->Get_ChiSq_Timing() + locChargedTrackHypothesis->Get_ChiSq_DCdEdx();
	double locFOM = (locNDF_Total > 0) ? TMath::Prob(locChiSq_Total, locNDF_Total) : numeric_limits<double>::quiet_NaN();
	locChargedTrackHypothesis->Set_ChiSq_Overall(locChiSq_Total, locNDF_Total, locFOM);

	locNewHypo->Set_ChiSq_Timing(double locChiSq, unsigned int locNDF);
	locNewHypo->Set_ChiSq_Overall(double locChiSq, unsigned int locNDF, double locFOM);

	return locNewHypo;
}


}

#endif // DParticleComboCreator_h
