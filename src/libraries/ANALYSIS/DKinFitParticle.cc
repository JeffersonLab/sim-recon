#include "DKinFitParticle.h"

DKinFitParticle::DKinFitParticle(void)
{
	Reset();
}

void DKinFitParticle::Reset(void)
{
	dPID = 0;
	dCharge = 0;
	dMass = 0.0;
	dSpacetimeVertex.SetXYZT(0.0, 0.0, 0.0, 0.0);
	dShowerEnergy = 0.0;
	dMomentum.SetXYZ(0.0, 0.0, 0.0);
	dCovarianceMatrix = NULL;
	dPathLength = 0.0;
	dPathLengthUncertainty = 0.0;

	dVertexConstraintFlag = 0;

	dKinFitParticleType = d_DetectedParticle;

	dEParamIndex = -1;
	dPxParamIndex = -1;
	dVxParamIndex = -1;
	dTParamIndex = -1;
	dLParamIndex = -1;

	dDecayingParticleAtProductionVertexFlag = true;

	dIsNeutralShowerFlag = false;
	dCommonVertexAndOrTimeConstraints.clear();
	dP4Constraints.clear();
}

const DKinFitConstraint_VertexBase* DKinFitParticle::Get_DefinedAtVertexAndOrTimeConstraint(void) const
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_VertexBase* locConstraint = dCommonVertexAndOrTimeConstraints[loc_i];
		deque<const DKinFitParticle*> locNoConstrainParticles = locConstraint->Get_NoConstrainParticles();
		for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
		{
			if(locNoConstrainParticles[loc_j] == this)
				return locConstraint;
		}
	}
	return NULL;
}

const DKinFitConstraint_VertexBase* DKinFitParticle::Get_ConstrainedAtVertexAndOrTimeConstraint(void) const
{
	for(size_t loc_i = 0; loc_i < dCommonVertexAndOrTimeConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_VertexBase* locConstraint = dCommonVertexAndOrTimeConstraints[loc_i];
		deque<const DKinFitParticle*> locFullConstrainParticles = locConstraint->Get_FullConstrainParticles();
		for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
		{
			if(locFullConstrainParticles[loc_j] == this)
				return locConstraint;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(locConstraint);
		if(locKinFitConstraint_Spacetime == NULL)
			continue;
		deque<const DKinFitParticle*> locOnlyConstrainTimeParticles = locKinFitConstraint_Spacetime->Get_OnlyConstrainTimeParticles();
		for(size_t loc_j = 0; loc_j < locOnlyConstrainTimeParticles.size(); ++loc_j)
		{
			if(locOnlyConstrainTimeParticles[loc_j] == this)
				return locConstraint;
		}
	}
	return NULL;
}

const DKinFitConstraint_P4* DKinFitParticle::Get_DefinedAtP4Constraint(void) const
{
	if(dP4Constraints.empty())
		return NULL;
	if((dKinFitParticleType == d_BeamParticle) || (dKinFitParticleType == d_DetectedParticle))
		return NULL;
	else if((dKinFitParticleType == d_TargetParticle) || (dKinFitParticleType == d_MissingParticle))
		return dP4Constraints[0];

	//Decaying Particle
	//true if the object's p3, v3, & t are defined at its production vertex (& common v & t are at decay vertex). else at it's decay vertex (& common v & t are at production vertex)
	bool locAtProductionVertexFlag = Get_DecayingParticleAtProductionVertexFlag();

	for(size_t loc_i = 0; loc_i < dP4Constraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locConstraint = dP4Constraints[loc_i];
		if(locAtProductionVertexFlag)
		{
			//search final state particles
			deque<const DKinFitParticle*> locParticles = locConstraint->Get_FinalParticles();
			for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
			{
				if(locParticles[loc_j] == this)
					return locConstraint;
			}
		}
		else
		{
			deque<const DKinFitParticle*> locParticles = locConstraint->Get_InitialParticles();
			for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
			{
				if(locParticles[loc_j] == this)
					return locConstraint;
			}
		}
	}
	return NULL;
}

const DKinFitConstraint_P4* DKinFitParticle::Get_ConstrainedAtP4Constraint(void) const
{
	if(dP4Constraints.empty())
		return NULL;
	if((dKinFitParticleType == d_BeamParticle) || (dKinFitParticleType == d_DetectedParticle))
		return dP4Constraints[0];
	else if((dKinFitParticleType == d_TargetParticle) || (dKinFitParticleType == d_MissingParticle))
		return NULL;

	//Decaying Particle
	//true if the object's p3, v3, & t are defined at its production vertex (& common v & t are at decay vertex). else at it's decay vertex (& common v & t are at production vertex)
	bool locAtProductionVertexFlag = Get_DecayingParticleAtProductionVertexFlag();

	for(size_t loc_i = 0; loc_i < dP4Constraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locConstraint = dP4Constraints[loc_i];
		if(locAtProductionVertexFlag)
		{
			//search initial state particles
			deque<const DKinFitParticle*> locParticles = locConstraint->Get_InitialParticles();
			for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
			{
				if(locParticles[loc_j] == this)
					return locConstraint;
			}
		}
		else
		{
			deque<const DKinFitParticle*> locParticles = locConstraint->Get_FinalParticles();
			for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
			{
				if(locParticles[loc_j] == this)
					return locConstraint;
			}
		}
	}
	return NULL;
}

const DKinFitConstraint_P4* DKinFitParticle::Get_P4ConstraintWhenInitial(void) const
{
	for(size_t loc_i = 0; loc_i < dP4Constraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locConstraint = dP4Constraints[loc_i];
		if((locConstraint->Get_InitialParticles())[0] == this)
			return locConstraint;
	}
	return NULL;
}

TVector3 DKinFitParticle::Get_CommonVertex(void) const
{
	const DKinFitConstraint_VertexBase* locConstraint = Get_ConstrainedAtVertexAndOrTimeConstraint();
	if(locConstraint != NULL)
		return locConstraint->Get_CommonVertex();
	locConstraint = Get_DefinedAtVertexAndOrTimeConstraint();
	return ((locConstraint != NULL) ? locConstraint->Get_CommonVertex() : TVector3());
}

double DKinFitParticle::Get_CommonTime(void) const
{
	const DKinFitConstraint_VertexBase* locConstraint = Get_ConstrainedAtVertexAndOrTimeConstraint();
	if(locConstraint != NULL)
	{
		const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
		return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_CommonTime() : 0.0);
	}
	locConstraint = Get_DefinedAtVertexAndOrTimeConstraint();
	if(locConstraint == NULL)
		return 0.0;
	const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
	return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_CommonTime() : 0.0);
}

TLorentzVector DKinFitParticle::Get_CommonSpacetimeVertex(void) const
{
	const DKinFitConstraint_VertexBase* locConstraint = Get_ConstrainedAtVertexAndOrTimeConstraint();
	if(locConstraint != NULL)
	{
		const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
		return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_CommonSpacetimeVertex() : (TLorentzVector(locConstraint->Get_CommonVertex(), 0.0)));
	}
	locConstraint = Get_DefinedAtVertexAndOrTimeConstraint();
	if(locConstraint == NULL)
		return (TLorentzVector());
	const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
	return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_CommonSpacetimeVertex() : (TLorentzVector(locConstraint->Get_CommonVertex(), 0.0)));
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

int DKinFitParticle::Get_CommonVxParamIndex(void) const
{
	const DKinFitConstraint_VertexBase* locConstraint = Get_ConstrainedAtVertexAndOrTimeConstraint();
	if(locConstraint != NULL)
		return locConstraint->Get_VxParamIndex();
	locConstraint = Get_DefinedAtVertexAndOrTimeConstraint();
	return ((locConstraint != NULL) ? locConstraint->Get_VxParamIndex() : -1);
}

int DKinFitParticle::Get_CommonTParamIndex(void) const
{
	const DKinFitConstraint_VertexBase* locConstraint = Get_ConstrainedAtVertexAndOrTimeConstraint();
	if(locConstraint != NULL)
	{
		const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
		return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_TParamIndex() : -1);
	}
	locConstraint = Get_DefinedAtVertexAndOrTimeConstraint();
	if(locConstraint == NULL)
		return -1;
	const DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<const DKinFitConstraint_Spacetime*>(locConstraint);
	return ((locSpacetimeConstraint != NULL) ? locSpacetimeConstraint->Get_TParamIndex() : -1);
}

