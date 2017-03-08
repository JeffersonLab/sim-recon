#include "ANALYSIS/DParticleCombo.h"


void DParticleCombo::Get_DecayChainParticles_Measured_Recursive(int locStepIndex, deque<const DKinematicData*>& locMeasuredParticles) const
{
	if((locStepIndex < 0) || (locStepIndex >= int(dParticleComboSteps.size())))
		return;
	const DParticleComboStep* locParticleComboStep = dParticleComboSteps[locStepIndex];

	if(locParticleComboStep->Get_InitialParticle_Measured() != NULL)
		locMeasuredParticles.push_back(locParticleComboStep->Get_InitialParticle_Measured());
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		if(locParticleComboStep->Is_FinalParticleDecaying(loc_i))
			Get_DecayChainParticles_Measured_Recursive(locParticleComboStep->Get_DecayStepIndex(loc_i), locMeasuredParticles);
		else if(locParticleComboStep->Is_FinalParticleDetected(loc_i))
			locMeasuredParticles.push_back(locParticleComboStep->Get_FinalParticle_Measured(loc_i));
	}
}

void DParticleCombo::Get_DetectedFinalChargedParticles_SourceObjects(deque<const DChargedTrack*>& locSourceChargedTracks) const
{
	locSourceChargedTracks.clear();
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboStep(loc_i)->Get_ParticleComboBlueprintStep();
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			const JObject* locJObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_j);
			if(locJObject == NULL)
				continue;
			const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locJObject);
			if(locChargedTrack == NULL)
				continue;
			locSourceChargedTracks.push_back(locChargedTrack);
		}
	}
}

void DParticleCombo::Get_DetectedFinalNeutralParticles_SourceObjects(deque<const DNeutralShower*>& locSourceNeutralShowers) const
{
	locSourceNeutralShowers.clear();
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboBlueprintStep* locParticleComboBlueprintStep = Get_ParticleComboStep(loc_i)->Get_ParticleComboBlueprintStep();
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_j)
		{
			const JObject* locJObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_j);
			if(locJObject == NULL)
				continue;
			const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locJObject);
			if(locNeutralShower == NULL)
				continue;
			locSourceNeutralShowers.push_back(locNeutralShower);
		}
	}
}

bool DParticleCombo::Check_AreMeasuredParticlesIdentical(const DParticleCombo* locParticleCombo) const
{
	if(locParticleCombo->Get_NumParticleComboSteps() != Get_NumParticleComboSteps())
		return false;

	deque<const DKinematicData*> locParticles, locCheckParticles;
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(loc_i);
		if(locParticleComboStep->Get_InitialParticle_Measured() != dParticleComboSteps[loc_i]->Get_InitialParticle_Measured())
			return false;
		if(locParticleComboStep->Get_TargetParticle() != dParticleComboSteps[loc_i]->Get_TargetParticle())
			return false;

		locParticleComboStep->Get_FinalParticles_Measured(locCheckParticles);
		dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locParticles);
		if(locParticles != locCheckParticles)
			return false;
	}
	return true;
}

set<pair<const JObject*, Particle_t> > DParticleCombo::Get_DecayingParticleSourceObjects(size_t locStepIndex) const
{
	set<pair<const JObject*, Particle_t> > locSourceObjects;
	const DParticleComboStep* locParticleComboStep = dParticleComboSteps[locStepIndex];
	for(size_t loc_i = 0; loc_i < locParticleComboStep->Get_NumFinalParticles(); ++loc_i)
	{
		//DecayStepIndex: one for each final particle: -2 if detected, -1 if missing, >= 0 if decaying, where the # is the step representing the particle decay
		int locDecayStepIndex = locParticleComboStep->Get_DecayStepIndex(loc_i);
		if(locDecayStepIndex == -2)
			locSourceObjects.insert(pair<const JObject*, Particle_t>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i), locParticleComboStep->Get_FinalParticleID(loc_i)));
		else if(locDecayStepIndex < 0)
			continue; //e.g. missing
		else
		{
			set<pair<const JObject*, Particle_t> > locNewSourceObjects = Get_DecayingParticleSourceObjects(locDecayStepIndex);
			locSourceObjects.insert(locNewSourceObjects.begin(), locNewSourceObjects.end());
		}
	}
	return locSourceObjects;
}

