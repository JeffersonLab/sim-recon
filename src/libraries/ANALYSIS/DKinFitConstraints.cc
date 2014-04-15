#include "DKinFitConstraints.h"

/******************************* DKinFitConstraint_VertexBase *******************************/

DKinFitConstraint_VertexBase::DKinFitConstraint_VertexBase(void)
{
	Reset();
}

void DKinFitConstraint_VertexBase::Reset(void)
{
	dVxParamIndex = -1;

	dDecayingParticles.clear();
	dFullConstrainParticles.clear();
	dNoConstrainParticles.clear();
	dDecayingParticlesToAssign.clear();
	dConstraintEquationParticleMap.clear();
}

int DKinFitConstraint_VertexBase::Get_FIndex(const DKinFitParticle* locKinFitParticle) const
{
	map<const DKinFitParticle*, unsigned int>::const_iterator locIterator = dConstraintEquationParticleMap.find(locKinFitParticle);
	if(locIterator == dConstraintEquationParticleMap.end())
		return -1;
	return locIterator->second;
}

void DKinFitConstraint_VertexBase::Replace_DecayingParticle(const DKinFitParticle* locOriginalParticle, const DKinFitParticle* locTreatAsDetectedParticle)
{
	for(deque<pair<DKinFitParticle*, bool> >::iterator locIterator = dDecayingParticles.begin(); locIterator != dDecayingParticles.end(); ++locIterator)
	{
		if((*locIterator).first != locOriginalParticle)
			continue;
		dDecayingParticles.erase(locIterator);
		break;
	}
	dDecayingParticlesToAssign.erase(const_cast<DKinFitParticle*>(locOriginalParticle));
	Add_FullConstrainParticle(const_cast<DKinFitParticle*>(locTreatAsDetectedParticle));
}

/******************************* DKinFitConstraint_Vertex *******************************/

DKinFitConstraint_Vertex::DKinFitConstraint_Vertex(void)
{
	Reset();
}

void DKinFitConstraint_Vertex::Reset(void)
{
	DKinFitConstraint_VertexBase::Reset();
	dCommonVertex.SetXYZ(0.0, 0.0, 0.0);
}

bool DKinFitConstraint_Vertex::Get_DecayingParticleInInitialStateFlag(const DKinFitParticle* locKinFitParticle) const
{
	for(size_t loc_i = 0; loc_i < dDecayingParticles.size(); ++loc_i)
	{
		if(dDecayingParticles[loc_i].first != locKinFitParticle)
			continue;
		return !(dDecayingParticles[loc_i].second);
	}
	return false;
}

/******************************* DKinFitConstraint_Spacetime *******************************/

DKinFitConstraint_Spacetime::DKinFitConstraint_Spacetime(void)
{
	Reset();
}

void DKinFitConstraint_Spacetime::Reset(void)
{
	DKinFitConstraint_VertexBase::Reset();

	dOnlyConstrainTimeParticles.clear();

	dUseRFTimeFlag = false;
	dCommonSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dTParamIndex = -1;
	dBeamParticle = NULL;
}

bool DKinFitConstraint_Spacetime::Get_DecayingParticleInInitialStateFlag(const DKinFitParticle* locKinFitParticle) const
{
	for(size_t loc_i = 0; loc_i < dDecayingParticles.size(); ++loc_i)
	{
		if(dDecayingParticles[loc_i].first != locKinFitParticle)
			continue;
		return !(dDecayingParticles[loc_i].second);
	}
	return false;
}

/******************************* DKinFitConstraint_P4 *******************************/

DKinFitConstraint_P4::DKinFitConstraint_P4(void)
{
	Reset();
}

void DKinFitConstraint_P4::Reset(void)
{
	dFIndex = 0;
	dConstrainedP4Particle = NULL;
	dInitialParticles.clear();
	dConstrainMassFlag = true;
	dConstrainMassByInvariantMassFlag = true;
	dIsActualP4ConstraintFlag = false;
	dConstrainedParticleIsInInitialStateFlag = true;
	dFinalParticles.clear();
}

void DKinFitConstraint_P4::Replace_Particle(const DKinFitParticle* locOriginalParticle, bool locInitialStateFlag, const DKinFitParticle* locTreatAsDetectedParticle)
{
	if(locInitialStateFlag)
	{
		for(size_t loc_i = 0; loc_i < dInitialParticles.size(); ++loc_i)
		{
			if(dInitialParticles[loc_i] != locOriginalParticle)
				continue;
			dInitialParticles[loc_i] = const_cast<DKinFitParticle*>(locTreatAsDetectedParticle);
			return;
		}
	}
	for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
	{
		if(dFinalParticles[loc_i] != locOriginalParticle)
			continue;
		dFinalParticles[loc_i] = const_cast<DKinFitParticle*>(locTreatAsDetectedParticle);
		return;
	}
}

