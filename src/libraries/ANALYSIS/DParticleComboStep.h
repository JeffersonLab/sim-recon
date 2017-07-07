#ifndef _DParticleComboStep_
#define _DParticleComboStep_

#include <vector>
#include <iostream>

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "ANALYSIS/DReactionStep.h"

using namespace std;
using namespace jana;

namespace DAnalysis
{

class DParticleComboStep
{
	public:
		// RESET:
		void Reset(void);

		//SET CONTENTS:
		void Set_Contents(const DKinematicData* locInitialParticle, const vector<const DKinematicData*>& locFinalParticles, const DLorentzVector& locSpacetimeVertex);

		// SET INITIAL PARTICLE:
		inline void Set_InitialParticle(const DKinematicData* locInitialParticle){dInitialParticle = locInitialParticle;}

		// SET FINAL PARTICLES:
		inline void Add_FinalParticle(const DKinematicData* locFinalParticle){dFinalParticles.push_back(locFinalParticle);}
		inline void Set_FinalParticle(const DKinematicData* locFinalParticle, size_t locFinalParticleIndex){dFinalParticles[locFinalParticleIndex] = locFinalParticle;}

		// SET MEASURED STEP:
		inline void Set_MeasuredParticleComboStep(const DParticleComboStep* locMeasuredParticleComboStep){dMeasuredStep = locMeasuredParticleComboStep;}

		// SET PRODUCTION/DECAY SPACETIME VERTEX
		inline void Set_SpacetimeVertex(const DLorentzVector& locSpacetimeVertex){dSpacetimeVertex = locSpacetimeVertex;}

		// GET INITIAL PARTICLES:
		inline const DKinematicData* Get_InitialParticle(void) const{return dInitialParticle;}
		const DKinematicData* Get_InitialParticle_Measured(void) const;

		// GET FINAL PARTICLES:
		inline size_t Get_NumFinalParticles(void) const{return dFinalParticles.size();}

		const DKinematicData* Get_FinalParticle(size_t locFinalParticleIndex) const;
		const DKinematicData* Get_FinalParticle_Measured(size_t locFinalParticleIndex) const;
		vector<const DKinematicData*> Get_FinalParticles(void) const{return dFinalParticles;}
		vector<const DKinematicData*> Get_FinalParticles_Measured(void) const{return ((dMeasuredStep == nullptr) ? dFinalParticles : dMeasuredStep->Get_FinalParticles());} //INCLUDES MISSING/DECAYING!!

		vector<const DKinematicData*> Get_FinalParticles(const DReactionStep* locReactionStep, bool locIncludeMissingFlag, bool locIncludeDecayingFlag = true, Charge_t locCharge = d_AllCharges) const;
		vector<const DKinematicData*> Get_FinalParticles_Measured(const DReactionStep* locReactionStep, Charge_t locCharge = d_AllCharges) const; //excludes missing/decaying!

		const JObject* Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const; //if missing or decaying: source object is nullptr!
		vector<const JObject*> Get_FinalParticle_SourceObjects(Charge_t locCharge = d_AllCharges) const;
		const DKinematicData* Get_MissingParticle(const DReactionStep* locReactionStep) const; //returns nullptr if none missing!

		// GET PRODUCTION/DECAY SPACETIME VERTEX
		inline DVector3 Get_Position(void) const{return dSpacetimeVertex.Vect();}
		inline double Get_Time(void) const{return dSpacetimeVertex.T();}
		inline DLorentzVector Get_SpacetimeVertex(void) const{return dSpacetimeVertex;}

	private:
		const DParticleComboStep* dMeasuredStep = nullptr;

		// INITIAL PARTICLES:
		const DKinematicData* dInitialParticle = nullptr; //if is nullptr: decaying or beam particle not yet set!

		// FINAL PARTICLES:
		vector<const DKinematicData*> dFinalParticles; //if particle is nullptr: missing or decaying! //these are DChargedTrackHypothesis or DNeutralParticleHypothesis objects if detected

		// PRODUCTION/DECAY SPACETIME VERTEX:
		DLorentzVector dSpacetimeVertex;
};

const JObject* Get_FinalParticle_SourceObject(const DKinematicData* locParticle);

inline const DKinematicData* DParticleComboStep::Get_InitialParticle_Measured(void) const
{
	return ((dMeasuredStep != nullptr) ? dMeasuredStep->Get_InitialParticle() : dInitialParticle);
}

inline void DParticleComboStep::Reset(void)
{
	dMeasuredStep = nullptr;
	dInitialParticle = nullptr;
	dFinalParticles.clear();
	dSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
}

inline void DParticleComboStep::Set_Contents(const DKinematicData* locInitialParticle, const vector<const DKinematicData*>& locFinalParticles, const DLorentzVector& locSpacetimeVertex)
{
	dInitialParticle = locInitialParticle;
	dFinalParticles = locFinalParticles;
	dSpacetimeVertex = locSpacetimeVertex;
}

inline const JObject* DParticleComboStep::Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const
{
	return DAnalysis::Get_FinalParticle_SourceObject(dFinalParticles[locFinalParticleIndex]);
}

inline const DKinematicData* DParticleComboStep::Get_FinalParticle(size_t locFinalParticleIndex) const
{
	if(locFinalParticleIndex >= dFinalParticles.size())
		return nullptr;
	return dFinalParticles[locFinalParticleIndex];
}

inline const DKinematicData* DParticleComboStep::Get_FinalParticle_Measured(size_t locFinalParticleIndex) const
{
	if(dMeasuredStep != nullptr)
		return dMeasuredStep->Get_FinalParticle_Measured(locFinalParticleIndex);
	if(locFinalParticleIndex >= dFinalParticles.size())
		return nullptr;
	return dFinalParticles[locFinalParticleIndex];
}

inline const DKinematicData* DParticleComboStep::Get_MissingParticle(const DReactionStep* locReactionStep) const
{
	int locMissingParticleIndex = locReactionStep->Get_MissingParticleIndex();
	if(locMissingParticleIndex == -1)
		return nullptr;
	return dFinalParticles[locMissingParticleIndex];
}

inline vector<const DKinematicData*> DParticleComboStep::Get_FinalParticles(const DReactionStep* locReactionStep, bool locIncludeMissingFlag, bool locIncludeDecayingFlag, Charge_t locCharge) const
{
	vector<const DKinematicData*> locFinalParticles;
	auto locMissingParticleIndex = locReactionStep->Get_MissingParticleIndex();
	for(int locPIDIndex = 0; locPIDIndex < int(locReactionStep->Get_NumFinalPIDs()); ++locPIDIndex)
	{
		if(!locIncludeMissingFlag && (locPIDIndex == locMissingParticleIndex))
			continue;
		if(!locIncludeDecayingFlag && (locPIDIndex != locMissingParticleIndex) && (Get_FinalParticle_SourceObject(locPIDIndex) == nullptr))
			continue;
		Particle_t locPID = locReactionStep->Get_FinalPID(locPIDIndex);
		if(Is_CorrectCharge(locPID, locCharge))
			locFinalParticles.push_back(dFinalParticles[locPIDIndex]);
	}
	return locFinalParticles;
}

inline vector<const DKinematicData*> DParticleComboStep::Get_FinalParticles_Measured(const DReactionStep* locReactionStep, Charge_t locCharge) const
{
	auto locStepPointer = (dMeasuredStep != nullptr) ? dMeasuredStep : this;
	return locStepPointer->Get_FinalParticles(locReactionStep, false, false, locCharge);
}

inline vector<const JObject*> DParticleComboStep::Get_FinalParticle_SourceObjects(Charge_t locCharge) const
{
	vector<const JObject*> locSourceObjects;
	for(size_t loc_i = 0; loc_i < Get_NumFinalParticles(); ++loc_i)
	{
		auto locSourceObject = Get_FinalParticle_SourceObject(loc_i);
		if(locSourceObject == nullptr)
			continue;
		if(Is_CorrectCharge(dFinalParticles[loc_i]->PID(), locCharge))
			locSourceObjects.push_back(locSourceObject);
	}
	return locSourceObjects;
}

inline vector<const DKinematicData*> Get_ParticlesWithPID(Particle_t locPID, const vector<const DKinematicData*>& locInputParticles)
{
	vector<const DKinematicData*> locOutputParticles;
	for(auto locParticle : locInputParticles)
	{
		if(locParticle == nullptr)
			continue;
		if(locParticle->PID() == locPID)
			locOutputParticles.push_back(locParticle);
	}
	return locOutputParticles;
}

inline const JObject* Get_FinalParticle_SourceObject(const DKinematicData* locParticle)
{
	if(locParticle == nullptr)
		return nullptr;

	auto locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(locParticle);
	if(locChargedHypo != nullptr)
	{
		const DChargedTrack* locChargedTrack = nullptr;
		locChargedHypo->GetSingle(locChargedTrack);
		return static_cast<const JObject*>(locChargedTrack);
	}

	auto locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(locParticle);
	return ((locNeutralHypo != nullptr) ? static_cast<const JObject*>(locNeutralHypo->Get_NeutralShower()) : nullptr);
}

} // end namespace

#endif // _DParticleComboStep_

