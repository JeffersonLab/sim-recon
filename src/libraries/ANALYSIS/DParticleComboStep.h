#ifndef _DParticleComboStep_
#define _DParticleComboStep_

#include <vector>
#include <iostream>

#include "particleType.h"
#include "DLorentzVector.h"
#include "PID/DKinematicData.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DChargedTrackHypothesis.h"

using namespace std;
using namespace jana;

class DParticleComboStep
{
	public:

		// CONSTRUCTOR:
		DParticleComboStep(void) : dMeasuredStep(nullptr), dInitialParticle(nullptr) {}

		// RESET:
		void Reset(void);

		//SET CONTENTS:
		void Set_Contents(const DKinematicData* locInitialParticle, const vector<const DKinematicData*>& locFinalParticles, const DLorentzVector& locSpacetimeVertex);;

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
		const JObject* Get_FinalParticle_SourceObject(size_t locFinalParticleIndex) const; //if missing or decaying: source object is nullptr!
		void Get_FinalParticles(deque<const DKinematicData*>& locParticles) const;
		const DKinematicData* Get_MissingParticle(void) const; //returns nullptr if none missing!

		// GET FINAL PARTICLES - BY PID:
		void Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const; //get all final particles of the given PID
		void Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const; //get all measured final particles of the given PID

		// GET FINAL PARTICLES - BY TRAIT:
		void Get_FinalParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const;
		void Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const;

		// GET PRODUCTION/DECAY SPACETIME VERTEX
		inline DVector3 Get_Position(void) const{return dSpacetimeVertex.Vect();}
		inline double Get_Time(void) const{return dSpacetimeVertex.T();}
		inline DLorentzVector Get_SpacetimeVertex(void) const{return dSpacetimeVertex;}

	private:
		const DParticleComboStep* dMeasuredStep;

		// INITIAL PARTICLES:
		const DKinematicData* dInitialParticle; //if is nullptr: decaying or beam particle not yet set!

		// FINAL PARTICLES:
		vector<const DKinematicData*> dFinalParticles; //if particle is nullptr: missing or decaying! //these are DChargedTrackHypothesis or DNeutralParticleHypothesis objects if detected

		// PRODUCTION/DECAY SPACETIME VERTEX:
		DLorentzVector dSpacetimeVertex;
};

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
	if(dFinalParticles[locFinalParticleIndex] == nullptr)
		return nullptr;

	auto locChargedHypo = dynamic_cast<const DChargedTrackHypothesis*>(dFinalParticles[locFinalParticleIndex]);
	if(locChargedHypo != nullptr)
		return static_cast<const JObject*>(locChargedHypo->Get_TrackTimeBased());

	auto locNeutralHypo = dynamic_cast<const DNeutralParticleHypothesis*>(dFinalParticles[locFinalParticleIndex]);
	return ((locNeutralHypo != nullptr) ? static_cast<const JObject*>(locNeutralHypo->Get_NeutralShower()) : nullptr);
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

inline void DParticleComboStep::Get_FinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		locParticles.push_back(dFinalParticles[loc_i]);
}

inline void DParticleComboStep::Get_FinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	if(dMeasuredStep != nullptr)
		dMeasuredStep->Get_FinalParticles_Measured(locParticles);
	else
	{
		for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

inline void DParticleComboStep::Get_FinalParticles(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(dFinalParticles[loc_i] == nullptr)
			continue;
		if(dFinalParticles[loc_i]->PID() == locPID)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

inline void DParticleComboStep::Get_FinalParticles_Measured(Particle_t locPID, deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	if(dMeasuredStep != nullptr)
		dMeasuredStep->Get_FinalParticles_Measured(locPID, locParticles);
	else
	{
		for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		{
			if(dFinalParticles[loc_i] == nullptr)
				continue;
			if(dFinalParticles[loc_i]->PID() == locPID)
				locParticles.push_back(dFinalParticles[loc_i]);
		}
	}
}

inline void DParticleComboStep::Get_DetectedFinalParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(Is_FinalParticleDetected(loc_i))
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

inline void DParticleComboStep::Get_DetectedFinalNeutralParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) == 0)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

inline void DParticleComboStep::Get_DetectedFinalChargedParticles(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(!Is_FinalParticleDetected(loc_i))
			continue;
		if(ParticleCharge(Get_FinalParticleID(loc_i)) != 0)
			locParticles.push_back(dFinalParticles[loc_i]);
	}
}

inline void DParticleComboStep::Get_DetectedFinalParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	if(dMeasuredStep != nullptr)
		dMeasuredStep->Get_DetectedFinalParticles_Measured(locParticles);
	else
	{
		for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		{
			if(Is_FinalParticleDetected(loc_i))
				locParticles.push_back(dFinalParticles[loc_i]);
		}
	}
}

inline void DParticleComboStep::Get_DetectedFinalNeutralParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	if(dMeasuredStep != nullptr)
		dMeasuredStep->Get_DetectedFinalNeutralParticles_Measured(locParticles);
	else
	{
		for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		{
			if(!Is_FinalParticleDetected(loc_i))
				continue;
			if(ParticleCharge(Get_FinalParticleID(loc_i)) == 0)
				locParticles.push_back(dFinalParticles[loc_i]);
		}
	}
}

inline void DParticleComboStep::Get_DetectedFinalChargedParticles_Measured(deque<const DKinematicData*>& locParticles) const
{
	locParticles.clear();
	if(dMeasuredStep != nullptr)
		dMeasuredStep->Get_DetectedFinalChargedParticles_Measured(locParticles);
	else
	{
		for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
		{
			if(!Is_FinalParticleDetected(loc_i))
				continue;
			if(ParticleCharge(Get_FinalParticleID(loc_i)) != 0)
				locParticles.push_back(dFinalParticles[loc_i]);
		}
	}
}

inline const DKinematicData* DParticleComboStep::Get_MissingParticle(void) const
{
	int locMissingParticleIndex = Get_MissingParticleIndex();
	if(locMissingParticleIndex == -1)
		return nullptr;
	return dFinalParticles[locMissingParticleIndex];
}


#endif // _DParticleComboStep_

