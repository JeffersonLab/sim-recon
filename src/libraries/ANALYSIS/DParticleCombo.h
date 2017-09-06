#ifndef _DParticleCombo_
#define _DParticleCombo_

#include <deque>
#include <vector>
#include <set>

#include "particleType.h"
#include "DResettable.h"
#include "PID/DKinematicData.h"
#include "PID/DEventRFBunch.h"
#include "ANALYSIS/DParticleComboStep.h"
#include "ANALYSIS/DKinFitResults.h"
#include "ANALYSIS/DReaction.h"

using namespace std;
using namespace DAnalysis;

namespace DAnalysis
{

class DParticleCombo : public DResettable
{
	public:
		void Reset(void);
		void Release(void){Reset();};

		// SET STEPS
		void Add_ParticleComboStep(const DParticleComboStep* locParticleComboStep){dParticleComboSteps.push_back(locParticleComboStep);}
		void Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex);

		// SET OBJECT DATA:
		void Set_KinFitResults(const DKinFitResults* locKinFitResults){dKinFitResults = locKinFitResults;}
		void Set_EventRFBunch(const DEventRFBunch* locEventRFBunch){dEventRFBunch = locEventRFBunch;}

		// GET OBJECT DATA:
		const DKinFitResults* Get_KinFitResults(void) const{return dKinFitResults;}
		const DEventRFBunch* Get_EventRFBunch(void) const{return dEventRFBunch;}

		// GET PARTCILE COMBO STEPS:
		size_t Get_NumParticleComboSteps(void) const{return dParticleComboSteps.size();}
		const DParticleComboStep* Get_ParticleComboStep(size_t locStepIndex) const;
		vector<const DParticleComboStep*> Get_ParticleComboSteps(void) const{return dParticleComboSteps;}

		// GET PARTICLES
		vector<const DKinematicData*> Get_MissingParticles(const DReaction* locReaction) const;
		vector<const DKinematicData*> Get_FinalParticles(const DReaction* locReaction, bool locIncludeMissingFlag = true, bool locIncludeDecayingFlag = true, Charge_t locCharge = d_AllCharges) const;
		vector<const DKinematicData*> Get_FinalParticles_Measured(const DReaction* locReaction, Charge_t locCharge = d_AllCharges) const;
		vector<const JObject*> Get_FinalParticle_SourceObjects(Charge_t locCharge = d_AllCharges) const;

		// GET FINAL PARTICLES - BY DECAY CHAIN:
		//get all of the measured particles included in the decaychain starting at locStepIndex
		vector<const DKinematicData*> Get_DecayChainParticles_Measured(const DReaction* locReaction, int locStepIndex) const;

		// OTHER:
		DLorentzVector Get_EventVertex(void) const;

	private:

		const DKinFitResults* dKinFitResults;
		const DEventRFBunch* dEventRFBunch;
		vector<const DParticleComboStep*> dParticleComboSteps;
};

inline void DParticleCombo::Reset(void)
{
	dKinFitResults = NULL;
	dEventRFBunch = NULL;
	dParticleComboSteps.clear();
}

inline const DParticleComboStep* DParticleCombo::Get_ParticleComboStep(size_t locStepIndex) const
{
	if(locStepIndex >= dParticleComboSteps.size())
		return NULL;
	return dParticleComboSteps[locStepIndex];
}

inline void DParticleCombo::Set_ParticleComboStep(const DParticleComboStep* locParticleComboStep, size_t locStepIndex)
{
	if(locStepIndex >= Get_NumParticleComboSteps())
		return;
	dParticleComboSteps[locStepIndex] = locParticleComboStep;
}

inline DLorentzVector DParticleCombo::Get_EventVertex(void) const
{
	if(dParticleComboSteps.empty())
		return DLorentzVector();
	return dParticleComboSteps[0]->Get_SpacetimeVertex();
}

inline vector<const DKinematicData*> DParticleCombo::Get_FinalParticles(const DReaction* locReaction, bool locIncludeMissingFlag, bool locIncludeDecayingFlag, Charge_t locCharge) const
{
	vector<const DKinematicData*> locParticles;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		auto locStepParticles = dParticleComboSteps[loc_i]->Get_FinalParticles(locReaction->Get_ReactionStep(loc_i), locIncludeMissingFlag, locIncludeDecayingFlag, locCharge);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
	return locParticles;
}

inline vector<const DKinematicData*> DParticleCombo::Get_FinalParticles_Measured(const DReaction* locReaction, Charge_t locCharge) const
{
	vector<const DKinematicData*> locParticles;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		auto locStepParticles = dParticleComboSteps[loc_i]->Get_FinalParticles_Measured(locReaction->Get_ReactionStep(loc_i), locCharge);
		locParticles.insert(locParticles.end(), locStepParticles.begin(), locStepParticles.end());
	}
	return locParticles;
}

inline vector<const JObject*> DParticleCombo::Get_FinalParticle_SourceObjects(Charge_t locCharge) const
{
	vector<const JObject*> locSourceObjects;
	for(size_t loc_i = 0; loc_i < Get_NumParticleComboSteps(); ++loc_i)
	{
		auto locStepSourceObjects = dParticleComboSteps[loc_i]->Get_FinalParticle_SourceObjects(locCharge);
		locSourceObjects.insert(locSourceObjects.end(), locStepSourceObjects.begin(), locStepSourceObjects.end());
	}
	return locSourceObjects;
}

inline vector<const DKinematicData*> DParticleCombo::Get_MissingParticles(const DReaction* locReaction) const
{
	vector<const DKinematicData*> locParticles;
	for(size_t loc_i = 0; loc_i < locReaction->Get_NumReactionSteps(); ++loc_i)
	{
		auto locMissingParticle = dParticleComboSteps[loc_i]->Get_MissingParticle(locReaction->Get_ReactionStep(loc_i));
		if(locMissingParticle != nullptr)
			locParticles.push_back(locMissingParticle);
	}
	return locParticles;
}

inline vector<const DKinematicData*> DParticleCombo::Get_DecayChainParticles_Measured(const DReaction* locReaction, int locStepIndex) const
{
	if((locStepIndex < 0) || (locStepIndex >= int(dParticleComboSteps.size())))
		return {};

	auto locParticleComboStep = dParticleComboSteps[locStepIndex];
	auto locReactionStep = locReaction->Get_ReactionStep(locStepIndex);
	auto locMissingIndex = locReactionStep->Get_MissingParticleIndex();

	vector<const DKinematicData*> locMeasuredParticles;
	if(locParticleComboStep->Get_InitialParticle_Measured() != NULL)
		locMeasuredParticles.push_back(locParticleComboStep->Get_InitialParticle_Measured());
	for(size_t loc_i = 0; loc_i < locReactionStep->Get_NumFinalPIDs(); ++loc_i)
	{
		if(int(loc_i) == locMissingIndex)
			continue;

		auto locMeasuredParticle = locParticleComboStep->Get_FinalParticle_Measured(loc_i);
		if(locMeasuredParticle == nullptr) //decaying!
		{
			auto locStepParticles = Get_DecayChainParticles_Measured(locReaction, DAnalysis::Get_DecayStepIndex(locReaction, locStepIndex, loc_i));
			locMeasuredParticles.insert(locMeasuredParticles.end(), locStepParticles.begin(), locStepParticles.end());
		}
		else
			locMeasuredParticles.push_back(locMeasuredParticle);
	}
	return locMeasuredParticles;
}

} // end namespace

#endif // _DParticleCombo_
