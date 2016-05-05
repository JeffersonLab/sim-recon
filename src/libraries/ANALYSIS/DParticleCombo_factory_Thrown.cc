// $Id$
//
//    File: DParticleCombo_factory_Thrown.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DParticleCombo_factory_Thrown.h"

//------------------
// init
//------------------
jerror_t DParticleCombo_factory_Thrown::init(void)
{
	MAX_dParticleComboStepPoolSize = 5;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DParticleCombo_factory_Thrown::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
 	vector<const DAnalysisUtilities*> locAnalysisUtilitiesVector;
	locEventLoop->Get(locAnalysisUtilitiesVector);
	dAnalysisUtilities = locAnalysisUtilitiesVector[0];
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DParticleCombo_factory_Thrown::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DParticleCombo_factory_Thrown::evnt()");
#endif

	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dParticleComboStepPool_All.size() > MAX_dParticleComboStepPoolSize)
	{
		for(size_t loc_i = MAX_dParticleComboStepPoolSize; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
			delete dParticleComboStepPool_All[loc_i];
		dParticleComboStepPool_All.resize(MAX_dParticleComboStepPoolSize);
	}
	dParticleComboStepPool_Available = dParticleComboStepPool_All;

	if(dParticleComboBlueprintStepPool_All.size() > MAX_dParticleComboStepPoolSize)
	{
		for(size_t loc_i = MAX_dParticleComboStepPoolSize; loc_i < dParticleComboBlueprintStepPool_All.size(); ++loc_i)
			delete dParticleComboBlueprintStepPool_All[loc_i];
		dParticleComboBlueprintStepPool_All.resize(MAX_dParticleComboStepPoolSize);
	}
	dParticleComboBlueprintStepPool_Available = dParticleComboBlueprintStepPool_All;

	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	dAnalysisUtilities->Get_ThrownParticleSteps(locEventLoop, locThrownSteps);

	if(locThrownSteps.empty())
		return NOERROR;

 	vector<const DReaction*> locReactions;
	locEventLoop->Get(locReactions, "Thrown");

	DParticleCombo* locParticleCombo = Build_ThrownCombo(locEventLoop, locReactions[0], locThrownSteps);
	_data.push_back(locParticleCombo);

	return NOERROR;
}

DParticleCombo* DParticleCombo_factory_Thrown::Build_ThrownCombo(JEventLoop* locEventLoop, const DReaction* locThrownReaction, deque<pair<const DMCThrown*, deque<const DMCThrown*> > >& locThrownSteps)
{
#ifdef VTRACE
	VT_TRACER("DParticleCombo_factory_Thrown::Build_ThrownCombo()");
#endif

 	vector<const DMCReaction*> locMCReactions;
	locEventLoop->Get(locMCReactions);

 	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Thrown");

	const DMCThrown* locMCThrown;
	DParticleCombo* locParticleCombo = new DParticleCombo();
	DParticleComboStep* locParticleComboStep = Get_ParticleComboStepResource();
	DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();
	locParticleCombo->Set_Reaction(locThrownReaction);
	locParticleCombo->Set_EventRFBunch(locEventRFBunches[0]);

	locParticleComboBlueprintStep->Set_ReactionStep(locThrownReaction->Get_ReactionStep(0));
	locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(-1);
	locParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboBlueprintStep);
	locParticleComboStep->Set_InitialParticle(&locMCReactions[0]->beam);
	locParticleComboStep->Set_TargetParticle(&locMCReactions[0]->target);
	locParticleComboStep->Set_Position(locMCReactions[0]->beam.position());
	locParticleComboStep->Set_Time(locMCReactions[0]->beam.time());

	for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
	{
		if(loc_i != 0) //else beam & target already set
		{
			locMCThrown = locThrownSteps[loc_i].first;
			locParticleComboStep = Get_ParticleComboStepResource();
			locParticleComboBlueprintStep = Get_ParticleComboBlueprintStepResource();

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
			locParticleComboBlueprintStep->Set_InitialParticleDecayFromStepIndex(locInitialParticleDecayFromStepIndex);

			locParticleComboBlueprintStep->Set_ReactionStep(locThrownReaction->Get_ReactionStep(loc_i));
			locParticleComboStep->Set_ParticleComboBlueprintStep(locParticleComboBlueprintStep);
			locParticleComboStep->Set_InitialParticle(locMCThrown);
			locParticleComboStep->Set_Position(locMCThrown->position());
			locParticleComboStep->Set_Time(locMCThrown->time());
		}
		for(size_t loc_j = 0; loc_j < locThrownSteps[loc_i].second.size(); ++loc_j)
		{
			locMCThrown = locThrownSteps[loc_i].second[loc_j];
			//check to see if this particle is a future decay parent
			int locDecayStepIndex = -2;
			for(size_t loc_k = loc_i + 1; loc_k < locThrownSteps.size(); ++loc_k)
			{
				if(locThrownSteps[loc_k].first != locMCThrown)
					continue;
				locDecayStepIndex = loc_k;
				break;
			}
			locParticleComboBlueprintStep->Add_FinalParticle_SourceObject(locMCThrown, locDecayStepIndex);
			locParticleComboStep->Add_FinalParticle(locMCThrown);
		}
		locParticleCombo->Add_ParticleComboStep(locParticleComboStep);
	}

	return locParticleCombo;
}

DParticleComboStep* DParticleCombo_factory_Thrown::Get_ParticleComboStepResource(void)
{
	DParticleComboStep* locParticleComboStep;
	if(dParticleComboStepPool_Available.empty())
	{
		locParticleComboStep = new DParticleComboStep();
		dParticleComboStepPool_All.push_back(locParticleComboStep);
	}
	else
	{
		locParticleComboStep = dParticleComboStepPool_Available.back();
		locParticleComboStep->Reset();
		dParticleComboStepPool_Available.pop_back();
	}
	return locParticleComboStep;
}

DParticleComboBlueprintStep* DParticleCombo_factory_Thrown::Get_ParticleComboBlueprintStepResource(void)
{
	DParticleComboBlueprintStep* locParticleComboBlueprintStep;
	if(dParticleComboBlueprintStepPool_Available.empty())
	{
		locParticleComboBlueprintStep = new DParticleComboBlueprintStep();
		dParticleComboBlueprintStepPool_All.push_back(locParticleComboBlueprintStep);
	}
	else
	{
		locParticleComboBlueprintStep = dParticleComboBlueprintStepPool_Available.back();
		locParticleComboBlueprintStep->Reset();
		dParticleComboBlueprintStepPool_Available.pop_back();
	}
	return locParticleComboBlueprintStep;
}

void DParticleCombo_factory_Thrown::Recycle_Combo(DParticleCombo* locParticleCombo)
{
	//deletes combo, but recycles steps
	for(size_t loc_i = 0; loc_i < locParticleCombo->Get_NumParticleComboSteps(); ++loc_i)
	{
		DParticleComboStep* locParticleComboStep = const_cast<DParticleComboStep*>(locParticleCombo->Get_ParticleComboStep(loc_i));
		DParticleComboBlueprintStep* locParticleComboBlueprintStep = const_cast<DParticleComboBlueprintStep*>(locParticleComboStep->Get_ParticleComboBlueprintStep());
		dParticleComboStepPool_Available.push_back(locParticleComboStep);
		dParticleComboBlueprintStepPool_Available.push_back(locParticleComboBlueprintStep);
	}
	delete locParticleCombo;
}

//------------------
// erun
//------------------
jerror_t DParticleCombo_factory_Thrown::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DParticleCombo_factory_Thrown::fini(void)
{
	for(size_t loc_i = 0; loc_i < dParticleComboStepPool_All.size(); ++loc_i)
		delete dParticleComboStepPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dParticleComboBlueprintStepPool_All.size(); ++loc_i)
		delete dParticleComboBlueprintStepPool_All[loc_i];

	return NOERROR;
}


