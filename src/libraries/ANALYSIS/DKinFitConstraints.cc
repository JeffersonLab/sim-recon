#include "DKinFitConstraints.h"

/******************************* DKinFitConstraint_VertexBase *******************************/

DKinFitConstraint_VertexBase::DKinFitConstraint_VertexBase(void)
{
	Reset();
}

void DKinFitConstraint_VertexBase::Reset(void)
{
	dVxParamIndex = -1;
}

int DKinFitConstraint_VertexBase::Get_FIndex(const DKinFitParticle* locKinFitParticle) const
{
	map<const DKinFitParticle*, unsigned int>::const_iterator locIterator = dConstraintEquationParticleMap.find(locKinFitParticle);
	if(locIterator == dConstraintEquationParticleMap.end())
		return -1;
	return locIterator->second;
}

/******************************* DKinFitConstraint_Vertex *******************************/

DKinFitConstraint_Vertex::DKinFitConstraint_Vertex(void)
{
	Reset();
}

void DKinFitConstraint_Vertex::Reset(void)
{
	DKinFitConstraint_VertexBase::Reset();

	dConstrainVertexParticles.clear();
	dDecayingParticles.clear();
	dNoConstrainParticles.clear();
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

	dConstrainSpacetimeParticles.clear();
	dOnlyConstrainTimeParticles.clear();
	dNoConstrainParticles.clear();
	dDecayingParticles.clear();

	dUseRFTimeFlag = false;
	dCommonSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dTParamIndex = -1;
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
	dFinalParticles.clear();
}

