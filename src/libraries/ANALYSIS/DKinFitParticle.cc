#include "DKinFitParticle.h"

DKinFitParticle::DKinFitParticle(void)
{
	Reset();
}

void DKinFitParticle::Reset(void)
{
	dCharge = 0;
	dMass = 0.0;
	dSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dShowerEnergy = 0.0;
	dMomentum.SetXYZ(0.0, 0.0, 0.0);
	dCovarianceMatrix = NULL;
	dPathLength = 0.0;

	dVertexConstraintFlag = 0;

	dKinFitParticleType = d_DetectedParticle;

	dEParamIndex = -1;
	dPxParamIndex = -1;
	dVxParamIndex = -1;
	dTParamIndex = -1;
	dLParamIndex = -1;

	dDecayingParticleAtProductionVertexFlag = true;

	dIsNeutralShowerFlag = false;
	dIsInP4FitFlag = false;
	dCommonVertexAndOrTimeConstraints.clear();
}

TVector3 DKinFitParticle::Get_CommonVertex(void) const
{
	DKinFitConstraint_VertexBase* locKinFitConstraint = Get_CommonVertexAndOrTimeConstraint();
	if(locKinFitConstraint == NULL)
		return (TVector3());
	return locKinFitConstraint->Get_CommonVertex();
}

double DKinFitParticle::Get_CommonTime(void) const
{
	DKinFitConstraint_VertexBase* locKinFitConstraint = Get_CommonVertexAndOrTimeConstraint();
	if(locKinFitConstraint == NULL)
		return 0.0;
	DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(locKinFitConstraint);
	if(locKinFitConstraint_Spacetime != NULL)
		return locKinFitConstraint_Spacetime->Get_CommonTime();
	return 0.0;
}

TLorentzVector DKinFitParticle::Get_CommonSpacetimeVertex(void) const
{
	DKinFitConstraint_VertexBase* locKinFitConstraint = Get_CommonVertexAndOrTimeConstraint();
	if(locKinFitConstraint == NULL)
		return (TLorentzVector());
	DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(locKinFitConstraint);
	if(locKinFitConstraint_Spacetime != NULL)
		return locKinFitConstraint_Spacetime->Get_CommonSpacetimeVertex();
	return (TLorentzVector());
}

bool DKinFitParticle::Get_IsInSpacetimeFitFlag(void) const
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
			return true;
	}
	return false;
}

bool DKinFitParticle::Get_IsInVertexFitFlag(void) const
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
			return true;
	}
	return false;
}

bool DKinFitParticle::Get_IsDefinedByVertexOrSpacetimeFitFlag(void) const //as opposed to helping constrain it OR not being in it
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			deque<DKinFitParticle*> locNoConstrainparticles = locKinFitConstraint_Spacetime->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainparticles.size(); ++loc_j)
			{
				if(locNoConstrainparticles[loc_j] == this)
					return true;
			}
			continue;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			deque<DKinFitParticle*> locNoConstrainparticles = locKinFitConstraint_Vertex->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainparticles.size(); ++loc_j)
			{
				if(locNoConstrainparticles[loc_j] == this)
					return true;
			}
			continue;
		}
	}
	return false;
}

bool DKinFitParticle::Get_IsConstrainedByVertexOrSpacetimeFitFlag(void) const //as opposed to being defined by it OR not being in it //INCLUDES neutral showers (constrained by time)
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			deque<DKinFitParticle*> locNoConstrainparticles = locKinFitConstraint_Spacetime->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainparticles.size(); ++loc_j)
			{
				if(locNoConstrainparticles[loc_j] == this)
					return false;
			}
			return true;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dCommonVertexAndOrTimeConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			deque<DKinFitParticle*> locNoConstrainparticles = locKinFitConstraint_Vertex->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainparticles.size(); ++loc_j)
			{
				if(locNoConstrainparticles[loc_j] == this)
					return false;
			}
			return true;
		}
	}
	return false;
}

int DKinFitParticle::Get_CommonVxParamIndex(void) const
{
	DKinFitConstraint_VertexBase* locKinFitConstraint = Get_CommonVertexAndOrTimeConstraint();
	return ((locKinFitConstraint == NULL) ? -1 : locKinFitConstraint->Get_VxParamIndex());
}

int DKinFitParticle::Get_CommonTParamIndex(void) const
{
	DKinFitConstraint_VertexBase* locKinFitConstraint = Get_CommonVertexAndOrTimeConstraint();
	if(locKinFitConstraint == NULL)
		return -1;
	DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(locKinFitConstraint);
	return ((locKinFitConstraint_Spacetime == NULL) ? -1 : locKinFitConstraint_Spacetime->Get_TParamIndex());
}

