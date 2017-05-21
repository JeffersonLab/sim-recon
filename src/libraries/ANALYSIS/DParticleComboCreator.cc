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


const DParticleCombo* DParticleComboCreator::Build_ParticleCombo(const DReactionVertexInfo* locReactionVertexInfo, const DSourceCombo* locFullCombo, const DKinematicData* locBeamParticle, int locRFBunchShift)
{
	auto locReaction = locReactionVertexInfo->Get_Reaction();
	auto locPreviousStepSourceCombo = locFullCombo;
	auto locParticleCombo = dResourcePool_ParticleCombo.Get_Resource();
	auto locIsPrimaryProductionVertex = locReactionVertexInfo->Get_StepVertexInfos().front()->Get_ProductionVertexFlag();

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
		DLorentzVector locSpacetimeVertex(locVertex, locTimeOffset + locEventRFBunch->dTime);

		//Get source objects, in order
		vector<const JObject*> locSourceObjects;
		auto locFinalPIDs = locReactionStep->Get_FinalPIDs();
		unordered_map<Particle_t, size_t> locPIDCountMap;
		for(size_t loc_j = 0; loc_j < locFinalPIDs.size(); ++loc_j)
		{
			if(loc_j == locReactionStep->Get_MissingParticleIndex())
			{
				locSourceObjects.push_back(nullptr);
				continue;
			}

			auto locPID = locFinalPIDs[loc_j];
			auto locPIDIterator = locPIDCountMap.find(locPID);
			if(locPIDIterator != locPIDCountMap.end())
				++(locPIDIterator->second);
			else
				locPIDCountMap.emplace(locPID, 1);

			size_t locPIDCountSoFar = 0;
			auto locSourceParticle = DAnalysis::Get_SourceParticle_ThisStep(locSourceCombo, locPID, locPIDCountMap[locPID], locPIDCountSoFar);
			locSourceObjects.push_back(locSourceParticle);
		}

		//Get/Build particles
		vector<const DKinematicData*> locFinalParticles;
		for(size_t loc_j = 0; loc_j < locSourceObjects.size(); ++loc_j)
		{
			auto locSourceObject = locSourceObjects[loc_j];
			if(locSourceObject == nullptr)
			{
				locFinalParticles.push_back(nullptr);
				continue;
			}

			auto locPID = locFinalPIDs[loc_j];
			if(ParticleCharge(locPID) == 0)
			{
				auto locNeutralShower = static_cast<const DNeutralShower*>(locNeutralShower);
			}
		}

		locParticleComboStep->Set_Contents(locReactionStep, locStepBeamParticle, locSourceObjects, locFinalParticles, locSpacetimeVertex);
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

const DNeutralParticleHypothesis* Create_NeutralHypothesis(JEventLoop* locEventLoop, const DNeutralShower* locNeutralShower, Particle_t locPID, const DEventRFBunch* locEventRFBunch, const DVector3 locVertex)
{
	//create
	DNeutralParticleHypothesis* locNeutralParticleHypothesis = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locEventLoop, locNeutralShower, locPID, locEventRFBunch, locVertex);
	if(locNeutralParticleHypothesis == NULL)
		continue;

	locNeutralParticleHypothesis->AddAssociatedObject(locNeutralParticles[loc_j]);

	const DNeutralParticleHypothesis* locOrigNeutralParticleHypothesis = locNeutralParticles[loc_j]->Get_Hypothesis(locPID);
	if(locOrigNeutralParticleHypothesis != NULL)
		locNeutralParticleHypothesis->AddAssociatedObject(locOrigNeutralParticleHypothesis);
}

void Create_Charged()
{
	//see if DChargedTrackHypothesis with the desired PID was created by the default factory, AND it passed the PreSelect cuts
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
	if(locChargedTrackHypothesis != NULL)
	{
		//it was/did. create new object with same PID (so that is registered with the combo factory, and because rf bunch could be different)
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
		DChargedTrackHypothesis* locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, dDetectorMatches, locEventRFBunch);

		locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
		locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
		locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrackHypothesis);
		_data.push_back(locNewChargedTrackHypothesis);
		return locNewChargedTrackHypothesis;
	}

}


}

#endif // DParticleComboCreator_h
