#include "ANALYSIS/DKinFitter_GlueX.h"

DKinFitter_GlueX::DKinFitter_GlueX(void)
{
	dIsBFieldNearBeamline = false;
	dMagneticFieldMap = NULL;
}

bool DKinFitter_GlueX::Get_IsBFieldNearBeamline(void) const
{
	return dIsBFieldNearBeamline;
}

void DKinFitter_GlueX::Set_BField(const DMagneticFieldMap* locMagneticFieldMap)
{
	dMagneticFieldMap = locMagneticFieldMap;
	dIsBFieldNearBeamline = (dMagneticFieldMap != NULL);
}

TVector3 DKinFitter_GlueX::Get_BField(const TVector3& locPosition) const
{
	if(dMagneticFieldMap == NULL)
		return (TVector3());

	double locBx, locBy, locBz;
	dMagneticFieldMap->GetField(locPosition.X(), locPosition.Y(),locPosition.Z(), locBx, locBy, locBz);
	return (TVector3(locBx, locBy, locBz));
}

void DKinFitter_GlueX::Reset_NewEvent(void)
{
	DKinFitter::Reset_NewEvent();
	dParticleMapping_InputToSource.clear();
}

void DKinFitter_GlueX::Reset_NewFit(void)
{
	DKinFitter::Reset_NewFit();
	dParticleMapping_OutputToSource.clear();
	dPulls.clear();
}

const DKinFitParticle* DKinFitter_GlueX::Make_BeamParticle(const DBeamPhoton* locBeamPhoton)
{
	const DKinematicData* locKinematicData = static_cast<const DKinematicData*>(locBeamPhoton);
	TLorentzVector locSpacetimeVertex(locKinematicData->position().X(),locKinematicData->position().Y(),locKinematicData->position().Z(), locKinematicData->time());
	TVector3 locMomentum(locKinematicData->momentum().X(),locKinematicData->momentum().Y(),locKinematicData->momentum().Z());
	Particle_t locPID = locKinematicData->PID();

	const DKinFitParticle* locKinFitParticle = DKinFitter::Make_BeamParticle(ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, &(locKinematicData->errorMatrix()));
	dParticleMapping_InputToSource[locKinFitParticle] = locKinematicData;
	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter_GlueX::Make_DetectedParticle(const DChargedTrackHypothesis* locChargedTrackHypothesis)
{
	const DKinematicData* locKinematicData = static_cast<const DKinematicData*>(locChargedTrackHypothesis);
	TLorentzVector locSpacetimeVertex(locKinematicData->position().X(),locKinematicData->position().Y(),locKinematicData->position().Z(), locKinematicData->time());
	TVector3 locMomentum(locKinematicData->momentum().X(),locKinematicData->momentum().Y(),locKinematicData->momentum().Z());
	Particle_t locPID = locKinematicData->PID();

	const DKinFitParticle* locKinFitParticle = DKinFitter::Make_DetectedParticle(ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, &(locKinematicData->errorMatrix()));
	dParticleMapping_InputToSource[locKinFitParticle] = locKinematicData;
	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter_GlueX::Make_DetectedParticle(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, bool locWillOnlyBeUsedInP4FitFlag)
{
	const DKinematicData* locKinematicData = static_cast<const DKinematicData*>(locNeutralParticleHypothesis);
	Particle_t locPID = locKinematicData->PID();

	const DKinFitParticle* locKinFitParticle;
	if(locWillOnlyBeUsedInP4FitFlag)
	{
		//use DNeutralParticleHypothesis object (assumes vertex is at target center! NOT IDEAL, AVOID IF POSSIBLE!!)
	  TLorentzVector locSpacetimeVertex(locKinematicData->position().X(),locKinematicData->position().Y(),locKinematicData->position().Z(), locKinematicData->time());
	  TVector3 locMomentum(locKinematicData->momentum().X(),locKinematicData->momentum().Y(),locKinematicData->momentum().Z());
	  locKinFitParticle = DKinFitter::Make_DetectedParticle(ParticleCharge(locPID), ParticleMass(locPID), locSpacetimeVertex, locMomentum, &(locKinematicData->errorMatrix()));
	}
	else
	{
		//use DNeutralShower object (doesn't make assumption about vertex!)
		vector<const DNeutralShower*> locNeutralShowers;
		locNeutralParticleHypothesis->GetT(locNeutralShowers);
		const DNeutralShower* locNeutralShower = locNeutralShowers[0];
		locKinFitParticle = DKinFitter::Make_DetectedShower(ParticleMass(locPID), TLorentzVector(locNeutralShower->dSpacetimeVertex.X(),locNeutralShower->dSpacetimeVertex.Y(),locNeutralShower->dSpacetimeVertex.Z(),locNeutralShower->dSpacetimeVertex.T()), locNeutralShower->dEnergy, &(locNeutralShower->dCovarianceMatrix));
	}
	dParticleMapping_InputToSource[locKinFitParticle] = locKinematicData;
	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter_GlueX::Make_DecayingParticle(Particle_t locPID)
{
	const DKinFitParticle* locKinFitParticle = DKinFitter::Make_DecayingParticle(ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMapping_InputToSource[locKinFitParticle] = NULL;
	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter_GlueX::Make_MissingParticle(Particle_t locPID)
{
	const DKinFitParticle* locKinFitParticle = DKinFitter::Make_MissingParticle(ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMapping_InputToSource[locKinFitParticle] = NULL;
	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter_GlueX::Make_TargetParticle(Particle_t locPID)
{
	const DKinFitParticle* locKinFitParticle = DKinFitter::Make_TargetParticle(ParticleCharge(locPID), ParticleMass(locPID));
	dParticleMapping_InputToSource[locKinFitParticle] = NULL;
	return locKinFitParticle;
}

bool DKinFitter_GlueX::Fit_Reaction(void)
{
	bool locFitStatus = DKinFitter::Fit_Reaction();
	if(!locFitStatus)
		return locFitStatus;

	const DKinFitParticle* locInputParticle;
	const DKinFitParticle* locOutputParticle;
	const DKinematicData* locKinematicData;

	dParticleMapping_OutputToSource.clear();
	dPulls.clear();

	//convert the map of the data: convert from (input particle -> data) to (output particle -> data)
	map<const DKinFitParticle*, const DKinFitParticle*> locKinFitParticleIOMap;
	DKinFitter::Get_KinFitParticleIOMap(locKinFitParticleIOMap); //map of input to output
	map<const DKinFitParticle*, const DKinematicData*>::const_iterator locParticleMappingIterator;
	for(locParticleMappingIterator = dParticleMapping_InputToSource.begin(); locParticleMappingIterator != dParticleMapping_InputToSource.end(); ++locParticleMappingIterator)
	{
		locInputParticle = locParticleMappingIterator->first;
		locKinematicData = locParticleMappingIterator->second;
		locOutputParticle = locKinFitParticleIOMap[locInputParticle];
		dParticleMapping_OutputToSource[locOutputParticle] = locKinematicData;
	}

	//convert the map of the pulls: convert from (output particle -> pull) to (kinematic data -> pull)
	map<const DKinFitParticle*, map<DKinFitPullType, double> > locPulls_KinFitter;
	DKinFitter::Get_Pulls(locPulls_KinFitter);
	map<const DKinFitParticle*, map<DKinFitPullType, double> >::iterator locPullIterator;
	for(locPullIterator = locPulls_KinFitter.begin(); locPullIterator != locPulls_KinFitter.end(); ++locPullIterator)
	{
		locOutputParticle = locPullIterator->first;
		locKinematicData = dParticleMapping_OutputToSource[locOutputParticle];
		dPulls[locKinematicData] = locPullIterator->second;
	}

	return locFitStatus;
}

bool DKinFitter_GlueX::Propagate_TrackInfoToCommonVertex(DKinematicData* locKinematicData, const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi)
{
	//locKinematicData must be generated from locKinFitParticle
		//this function should only be used on decaying particles involved in two vertex fits:
			//propagates the track information from the vertex at which it is DEFINED to the OTHER vertex (e.g. production -> decay)
	TVector3 locMomentum;
	TLorentzVector locSpacetimeVertex;
	pair<double, double> locPathLengthPair;
	TMatrixDSym* locCovarianceMatrix = DKinFitter::Get_MatrixDSymResource();
	if(!DKinFitter::Propagate_TrackInfoToCommonVertex(locKinFitParticle, locVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, *locCovarianceMatrix))
		return false;
	locKinematicData->setMomentum(DVector3(locMomentum.X(),locMomentum.Y(),locMomentum.Z()));
	locKinematicData->setPosition(DVector3(locSpacetimeVertex.Vect().X(),locSpacetimeVertex.Vect().Y(),locSpacetimeVertex.Vect().Z()));
	locKinematicData->setTime(locSpacetimeVertex.T());
	locKinematicData->setErrorMatrix(*locCovarianceMatrix);
	locKinematicData->setPathLength(locPathLengthPair.first, locPathLengthPair.second);
	return true;
}



