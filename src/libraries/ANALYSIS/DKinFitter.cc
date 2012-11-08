#include "DKinFitter.h"

//THINGS TO DO:
	//double check: spacetime eqs & partial derivatives
	//double check: Spacetime being set for noconstrain particles
	//double check: vertex constraints for linked decaying particles being set correctly
	//double check: propagation of output covariance matrix, track time, and path lengths done correctly
	//partial derivatives
		//think about allowing missing beam particle in rf time constraint
		//and think about correlations between beam particle and rf time

//CONSTRAINTS:
	//Constraint types: p4 conservation, common vertex, common spacetime
	//P4 Conservation Notes:
		//P4 is required to be conserved between the particles at the common vertex (if it is simultaneously kinematically fit)
		//For input neutral showers, the momentum is defined-by and updated-with the common vertex (see "Common Vertex Notes")
			//if a common vertex is not simultaneously kinematically fit for a neutral shower, the fit will return false
			//To include a neutral shower in a p4 fit without a vertex fit, create a particle from it with a full 7x7 covaraince matrix and input it instead of the neutral shower
		//If the particle is charged, the magnetic field is taken into account if a common vertex is simultaneously kinematically fit (magnetic field assumed to be constant over the propagation range)
			//If a common vertex is not defined, the momentum is assumed to be constant with position
		//Energy loss between the common vertex and the particle position is ignored
	//Common Vertex Notes:
		//See "P4 Conservation Notes" section about vertices being used in those constraints.
		//If the particle is charged, the magnetic field is taken into account (magnetic field assumed to be constant over the propagation range)
		//You can include the beam particle in the common vertex fit.  However, it is ignored if it's xy uncertainties and/or(?) xy parameters are zero (matrix won't invert).  
		//Neutral and missing particles included in the constraint will not be used to constrain the vertex, but will be set with the fit common vertex
			//This is necessary if you want the neutral particle momentum (from an input shower) to change with the reconstructed vertex
		//If a decaying particle is set in two vertex constraints (even if one is a spacetime constraint):
			//It's position will be defined by the point at which it's p4 is constrained, and it will be constrained to have a common vertex at the other vertex
	//Common Spacetime Notes:
		//It is not possible to fit a common time without simultaneously fitting a common vertex.
			//Requiring a common time at a non-common vertex has no meaning, and the fitter is not setup to input a pre-fit common vertex with uncertainties.
		//If the particle is charged, the magnetic field is taken into account (magnetic field assumed to be constant over the propagation range)
		//You can include both the beam particle AND the RF time in the common time fit
			//the RF parameters are defined to be the same as those of the beam particle, except the time and it's uncertainty
		//Missing particles included in the constraint will not be used to constrain the spacetime, but will be set with the fit common spacetime
		//If a decaying particle is set in two spacetime constraints:
			//It's spacetime will be defined by the point at which it's p4 is constrained, and it will be constrained to have a common spacetime at the other vertex

//SETUP:
	//If setup of a particle or constraint fails, it will return NULL instead of the created object
	//only input a neutral shower if its vertex will be defined by a simultaneous vertex fit!
		//else import it as a full detected particle with a 7x7 covariance matrix

DKinFitConstraint::~DKinFitConstraint(void){}

DKinFitter::DKinFitter(void)
{
	dLinkVerticesFlag = false;
	dDebugLevel = 0;

	dMaxKinFitParticlePoolSize = 100;
	dMaxKinFitConstraintVertexPoolSize = 25;
	dMaxKinFitConstraintSpacetimePoolSize = 25;
	dMaxKinFitConstraintP4PoolSize = 25;
	dMaxMatrixDSymPoolSize = 100;

	dMaxNumIterations = 10;

	Reset_NewEvent();
}

DKinFitter::~DKinFitter(void)
{
	for(size_t loc_i = 0; loc_i < dKinFitParticlePool_All.size(); ++loc_i)
		delete dKinFitParticlePool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintVertexPool_All.size(); ++loc_i)
		delete dKinFitConstraintVertexPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintSpacetimePool_All.size(); ++loc_i)
		delete dKinFitConstraintSpacetimePool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintP4Pool_All.size(); ++loc_i)
		delete dKinFitConstraintP4Pool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dMatrixDSymPool_All.size(); ++loc_i)
		delete dMatrixDSymPool_All[loc_i];
}

void DKinFitter::Reset_NewEvent(void)
{
	Reset_NewFit();

	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dKinFitParticlePool_All.size() > dMaxKinFitParticlePoolSize)
	{
		for(size_t loc_i = dMaxKinFitParticlePoolSize; loc_i < dKinFitParticlePool_All.size(); ++loc_i)
			delete dKinFitParticlePool_All[loc_i];
		dKinFitParticlePool_All.resize(dMaxKinFitParticlePoolSize);
	}
	dKinFitParticlePool_Available = dKinFitParticlePool_All;

	if(dKinFitConstraintVertexPool_All.size() > dMaxKinFitConstraintVertexPoolSize)
	{
		for(size_t loc_i = dMaxKinFitConstraintVertexPoolSize; loc_i < dKinFitConstraintVertexPool_All.size(); ++loc_i)
			delete dKinFitConstraintVertexPool_All[loc_i];
		dKinFitConstraintVertexPool_All.resize(dMaxKinFitConstraintVertexPoolSize);
	}
	dKinFitConstraintVertexPool_Available = dKinFitConstraintVertexPool_All;

	if(dKinFitConstraintSpacetimePool_All.size() > dMaxKinFitConstraintSpacetimePoolSize)
	{
		for(size_t loc_i = dMaxKinFitConstraintSpacetimePoolSize; loc_i < dKinFitConstraintSpacetimePool_All.size(); ++loc_i)
			delete dKinFitConstraintSpacetimePool_All[loc_i];
		dKinFitConstraintSpacetimePool_All.resize(dMaxKinFitConstraintSpacetimePoolSize);
	}
	dKinFitConstraintSpacetimePool_Available = dKinFitConstraintSpacetimePool_All;

	if(dKinFitConstraintP4Pool_All.size() > dMaxKinFitConstraintP4PoolSize)
	{
		for(size_t loc_i = dMaxKinFitConstraintP4PoolSize; loc_i < dKinFitConstraintP4Pool_All.size(); ++loc_i)
			delete dKinFitConstraintP4Pool_All[loc_i];
		dKinFitConstraintP4Pool_All.resize(dMaxKinFitConstraintP4PoolSize);
	}
	dKinFitConstraintP4Pool_Available = dKinFitConstraintP4Pool_All;

	if(dMatrixDSymPool_All.size() > dMaxMatrixDSymPoolSize)
	{
		for(size_t loc_i = dMaxMatrixDSymPoolSize; loc_i < dMatrixDSymPool_All.size(); ++loc_i)
			delete dMatrixDSymPool_All[loc_i];
		dMatrixDSymPool_All.resize(dMaxMatrixDSymPoolSize);
	}
	dMatrixDSymPool_Available = dMatrixDSymPool_All;
}

void DKinFitter::Reset_NewFit(void)
{
	dKinFitConstraints.clear();
	dKinFitParticles.clear();
	dKinFitParticleIOMap.clear();

	dRFMatchedBeamParticle = NULL;
	dRFTimeParamIndex = -1;
	dRFTime = 0.0;
	dRFUncertainty = 0.0;

	dNumXi = 0;
	dNumEta = 0;
	dNumF = 0;

	dChiSq = 0.0;
	dNDF = 0;
	dConfidenceLevel = 0.0;
	dPulls.clear();

	dVXi = Get_MatrixDSymResource();
	dVEta = Get_MatrixDSymResource();
}

DKinFitParticle* DKinFitter::Get_KinFitParticleResource(void)
{
	DKinFitParticle* locKinFitParticle;
	if(dKinFitParticlePool_Available.empty())
	{
		locKinFitParticle = new DKinFitParticle;
		dKinFitParticlePool_All.push_back(locKinFitParticle);
	}
	else
	{
		locKinFitParticle = dKinFitParticlePool_Available.back();
		locKinFitParticle->Reset();
		dKinFitParticlePool_Available.pop_back();
	}
	return locKinFitParticle;
}

DKinFitConstraint_Vertex* DKinFitter::Get_KinFitConstraintVertexResource(void)
{
	DKinFitConstraint_Vertex* locKinFitConstraint;
	if(dKinFitConstraintVertexPool_Available.empty())
	{
		locKinFitConstraint = new DKinFitConstraint_Vertex;
		dKinFitConstraintVertexPool_All.push_back(locKinFitConstraint);
	}
	else
	{
		locKinFitConstraint = dKinFitConstraintVertexPool_Available.back();
		locKinFitConstraint->Reset();
		dKinFitConstraintVertexPool_Available.pop_back();
	}
	return locKinFitConstraint;
}

DKinFitConstraint_Spacetime* DKinFitter::Get_KinFitConstraintSpacetimeResource(void)
{
	DKinFitConstraint_Spacetime* locKinFitConstraint;
	if(dKinFitConstraintSpacetimePool_Available.empty())
	{
		locKinFitConstraint = new DKinFitConstraint_Spacetime;
		dKinFitConstraintSpacetimePool_All.push_back(locKinFitConstraint);
	}
	else
	{
		locKinFitConstraint = dKinFitConstraintSpacetimePool_Available.back();
		locKinFitConstraint->Reset();
		dKinFitConstraintSpacetimePool_Available.pop_back();
	}
	return locKinFitConstraint;
}

DKinFitConstraint_P4* DKinFitter::Get_KinFitConstraintP4Resource(void)
{
	DKinFitConstraint_P4* locKinFitConstraint;
	if(dKinFitConstraintP4Pool_Available.empty())
	{
		locKinFitConstraint = new DKinFitConstraint_P4;
		dKinFitConstraintP4Pool_All.push_back(locKinFitConstraint);
	}
	else
	{
		locKinFitConstraint = dKinFitConstraintP4Pool_Available.back();
		locKinFitConstraint->Reset();
		dKinFitConstraintP4Pool_Available.pop_back();
	}
	return locKinFitConstraint;
}

TMatrixDSym* DKinFitter::Get_MatrixDSymResource(void)
{
	TMatrixDSym* locMatrixDSym;
	if(dMatrixDSymPool_Available.empty())
	{
		locMatrixDSym = new TMatrixDSym();
		dMatrixDSymPool_All.push_back(locMatrixDSym);
	}
	else
	{
		locMatrixDSym = dMatrixDSymPool_Available.back();
		dMatrixDSymPool_Available.pop_back();
	}
	return locMatrixDSym;
}

TMatrixDSym* DKinFitter::Clone_MatrixDSym(const TMatrixDSym* locMatrix)
{
	if(locMatrix == NULL)
		return NULL;
	TMatrixDSym* locNewMatrix = Get_MatrixDSymResource();
	int locMatrixSize = locMatrix->GetNcols();
	locNewMatrix->ResizeTo(locMatrixSize, locMatrixSize);
	*locNewMatrix = *locMatrix;
	return locNewMatrix;
}

void DKinFitter::Clone_KinFitParticles(const deque<const DKinFitParticle*>& locKinFitParticles, deque<DKinFitParticle*>& locClonedKinFitParticles)
{
	locClonedKinFitParticles.clear();
	for(size_t loc_i = 0; loc_i < locKinFitParticles.size(); ++loc_i)
		locClonedKinFitParticles.push_back(Clone_KinFitParticle(locKinFitParticles[loc_i]));
}

DKinFitParticle* DKinFitter::Clone_KinFitParticle(const DKinFitParticle* locKinFitParticle)
{
	DKinFitParticle* locClonedKinFitParticle = Get_KinFitParticleResource();
	*locClonedKinFitParticle = *locKinFitParticle;

	const TMatrixDSym* locCovarianceMatrix = locKinFitParticle->Get_CovarianceMatrix();
	if(locCovarianceMatrix != NULL)
		locClonedKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	return locClonedKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_DecayingParticle(int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_DecayingParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Decaying particle set. Q, Mass = " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_MissingParticle(int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_MissingParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Missing particle set. Q, Mass = " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_BeamParticle(int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_BeamParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Beam particle set. Q, Mass, P3, V3, T = " << locCharge << ", " << locMass << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_TargetParticle(int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_TargetParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Target particle set. Q, Mass = " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_DetectedParticle(int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Detected particle set. Q, Mass, P3, V3, T = " << locCharge << ", " << locMass << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_DetectedShower(double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 5) || (locCovarianceMatrix->GetNcols() != 5))
		return NULL; //is not 5x5

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_Charge(0);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_IsNeutralShowerFlag(true);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_ShowerEnergy(locShowerEnergy);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Detected shower set. Q, Mass, E, V3, T = 0, " << locMass << ", " << locShowerEnergy << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

	return locKinFitParticle;
}

void DKinFitter::Set_RFTime(double locRFTime, double locRFUncertainty, const DKinFitParticle* locRFMatchedBeamParticle)
{
	cout << "ERROR: THIS IS NOT SUPPORTED YET. RETURNING." << endl;
	return;

	dRFTime = locRFTime;
	dRFUncertainty = locRFUncertainty;
	dRFMatchedBeamParticle = GetOrCreate_ClonedParticle(locRFMatchedBeamParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: RF Time set. t = " << locRFTime << endl;
}

DKinFitParticle* DKinFitter::GetOrCreate_ClonedParticle(const DKinFitParticle* locKinFitParticle)
{
	map<const DKinFitParticle*, DKinFitParticle*>::iterator locIterator = dKinFitParticleIOMap.find(locKinFitParticle);
	if(locIterator != dKinFitParticleIOMap.end())
		return locIterator->second;
	DKinFitParticle* locClonedKinFitParticle = Clone_KinFitParticle(locKinFitParticle);
	dKinFitParticleIOMap[locKinFitParticle] = locClonedKinFitParticle;
	dKinFitParticles.push_back(locClonedKinFitParticle);
	return locClonedKinFitParticle;
}

const DKinFitConstraint_Vertex* DKinFitter::Add_VertexConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, TVector3& locVertexGuess)
{
	deque<DKinFitParticle*> locConstrainVertexParticles; //charged particles, decaying particles, beam particles
	deque<pair<DKinFitParticle*, bool> > locDecayingParticles; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
	deque<DKinFitParticle*> locNoConstrainParticles; //missing particles & neutral showers //not used to constrain vertex or time, but fit vertex is set for this particle

	//require all of the tracks to pass through a common point
	//decaying particles are only used to constrain the fit if they are included in exactly two vertex constraints (once as an initial particle, once as a final particle)
		//else they are treated as dNoConstrainParticles
	//will set decaying particles in dConstrainVertexParticles and dNoConstrainParticles when ready, but not yet!

	//first check to make sure the inputs are ok //can't tell if enough particles until Resolve_DecayingParticleSpacetimeLinks() //due to decaying particles in 2 constraints
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		if(locInitialParticles[loc_i] == NULL)
			return NULL;
		else if(locInitialParticles[loc_i]->Get_KinFitParticleType() == d_DetectedParticle)
			return NULL;
	}
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		if(locFinalParticles[loc_i] == NULL)
			return NULL;
		else if(locFinalParticles[loc_i]->Get_KinFitParticleType() == d_BeamParticle)
			return NULL;
		else if(locFinalParticles[loc_i]->Get_KinFitParticleType() == d_TargetParticle)
			return NULL;
	}

	//constraint will be created.  first clone the particles (if not previously cloned)
	deque<DKinFitParticle*> locClonedInitialParticles, locClonedFinalParticles;
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
		locClonedInitialParticles.push_back(GetOrCreate_ClonedParticle(locInitialParticles[loc_i]));
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
		locClonedFinalParticles.push_back(GetOrCreate_ClonedParticle(locFinalParticles[loc_i]));

	//sort particles by how they'll be used by the constraint
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	{
		if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_TargetParticle)
			locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		else if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locClonedInitialParticles[loc_i], false));
			else
				locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		}
		else if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_BeamParticle)
		{
			//only add if both vx & vy uncertainties are non-zero (else constraints are bad!!)
			const TMatrixDSym& locCovarianceMatrix = *(locClonedInitialParticles[loc_i]->Get_CovarianceMatrix());
			if((locCovarianceMatrix(3, 3) > 0.0) && (locCovarianceMatrix(4, 4) > 0.0)) //include beamline in vertex fit!
				locConstrainVertexParticles.push_back(locClonedInitialParticles[loc_i]);
			else
				locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		}
	}
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	{
		if(locClonedFinalParticles[loc_i]->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locClonedFinalParticles[loc_i]);
		else if(locClonedFinalParticles[loc_i]->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locClonedFinalParticles[loc_i], true));
			else
				locNoConstrainParticles.push_back(locClonedFinalParticles[loc_i]);
		}
		else if(locClonedFinalParticles[loc_i]->Get_Charge() == 0)
			locNoConstrainParticles.push_back(locClonedFinalParticles[loc_i]);
		else
			locConstrainVertexParticles.push_back(locClonedFinalParticles[loc_i]);
	}

	//set vertex constraint flags
	TVector3 locMomentum;
	unsigned short int locVertexConstraintFlag;
	for(size_t loc_i = 0; loc_i < locConstrainVertexParticles.size(); ++loc_i)
	{
		locMomentum = locConstrainVertexParticles[loc_i]->Get_Momentum();
		if(fabs(locMomentum.Pz()) > fabs(locMomentum.Px()))
			locVertexConstraintFlag = (fabs(locMomentum.Pz()) > fabs(locMomentum.Py())) ? 1 : 2;
		else
			locVertexConstraintFlag = (fabs(locMomentum.Px()) > fabs(locMomentum.Py())) ? 3 : 2;
		locConstrainVertexParticles[loc_i]->Set_VertexConstraintFlag(locVertexConstraintFlag);
	}

	//set momentum of neutral showers that have 5x5 covariance matrix (shower energy input instead of p3)
	for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	{
		if(!locNoConstrainParticles[loc_i]->Get_IsNeutralShowerFlag())
			continue; //only do for neutral showers

		double locE = locNoConstrainParticles[loc_i]->Get_ShowerEnergy();
		double locMass = locNoConstrainParticles[loc_i]->Get_Mass();
		double locPMag = sqrt(locE*locE - locMass*locMass);
		TVector3 locMomentum = locNoConstrainParticles[loc_i]->Get_Position() - locNoConstrainParticles[loc_i]->Get_CommonVertex();
		locMomentum.SetMag(locPMag);
		locNoConstrainParticles[loc_i]->Set_Momentum(locMomentum);
	}

	//create the constraint and set its members
	DKinFitConstraint_Vertex* locKinFitConstraint = Get_KinFitConstraintVertexResource();
	locKinFitConstraint->Set_ConstrainVertexParticles(locConstrainVertexParticles);
	locKinFitConstraint->Set_DecayingParticles(locDecayingParticles);
	locKinFitConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locKinFitConstraint->Set_CommonVertex(locVertexGuess);

	//add constraint to particles
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
		locClonedInitialParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locKinFitConstraint);
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
		locClonedFinalParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locKinFitConstraint);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Vertex constraint added. Constrained particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locConstrainVertexParticles.size(); ++loc_i)
	 		cout << locConstrainVertexParticles[loc_i]->Get_Charge() << ", " << locConstrainVertexParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	 		cout << locNoConstrainParticles[loc_i]->Get_Charge() << ", " << locNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	 		cout << locDecayingParticles[loc_i].first->Get_Charge() << ", " << locDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locKinFitConstraint));
	return locKinFitConstraint;
}

const DKinFitConstraint_Spacetime* DKinFitter::Add_SpacetimeConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locUseRFTimeFlag, TVector3& locVertexGuess, double locCommonTimeGuess)
{
	cout << "ERROR: THIS IS NOT SUPPORTED YET. RETURNING." << endl;
	return NULL;

	deque<DKinFitParticle*> locConstrainSpacetimeParticles; //charged particles, decaying particles, beam particles
	deque<DKinFitParticle*> locOnlyConstrainTimeParticles; //neutral showers //not used to constrain vertex, but fit vertex is used for time constraint
	deque<pair<DKinFitParticle*, bool> > locDecayingParticles; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
	deque<DKinFitParticle*> locNoConstrainParticles; //missing particles //not used to constrain vertex or time, but fit vertex & time are set for this particle

	//require all of the tracks to pass through a common point at a common time
	//decaying particles are only used to constrain the fit if they are included in exactly two vertex constraints (once as an initial particle, once as a final particle)
		//else they are treated as dNoConstrainParticles
	//will set decaying particles in dConstrainSpacetimeParticles and dNoConstrainParticles when ready, but not yet!

	//first check to make sure the inputs are ok //can't tell if enough particles until Resolve_DecayingParticleSpacetimeLinks() //due to decaying particles in 2 constraints
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		if(locInitialParticles[loc_i] == NULL)
			return NULL;
		else if(locInitialParticles[loc_i]->Get_KinFitParticleType() == d_DetectedParticle)
			return NULL;
	}
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		if(locFinalParticles[loc_i] == NULL)
			return NULL;
		else if(locFinalParticles[loc_i]->Get_KinFitParticleType() == d_BeamParticle)
			return NULL;
		else if(locFinalParticles[loc_i]->Get_KinFitParticleType() == d_TargetParticle)
			return NULL;
	}

	//constraint will be created.  first clone the particles (if not previously cloned)
	deque<DKinFitParticle*> locClonedInitialParticles, locClonedFinalParticles;
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
		locClonedInitialParticles.push_back(GetOrCreate_ClonedParticle(locInitialParticles[loc_i]));
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
		locClonedFinalParticles.push_back(GetOrCreate_ClonedParticle(locFinalParticles[loc_i]));

	//sort particles by how they'll be used by the constraint
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	{
		if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		else if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_TargetParticle)
			locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		else if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locClonedInitialParticles[loc_i], false));
			else
				locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		}
		else if(locClonedInitialParticles[loc_i]->Get_KinFitParticleType() == d_BeamParticle)
		{
			//only add if both vx & vy uncertainties are non-zero (else constraints are bad!!)
			const TMatrixDSym& locCovarianceMatrix = *(locClonedInitialParticles[loc_i]->Get_CovarianceMatrix());
			if((locCovarianceMatrix(3, 3) > 0.0) && (locCovarianceMatrix(4, 4) > 0.0))
				locConstrainSpacetimeParticles.push_back(locClonedInitialParticles[loc_i]); //include beamline in vertex fit!
			else
				locNoConstrainParticles.push_back(locClonedInitialParticles[loc_i]);
		}
	}
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	{
		if(locClonedFinalParticles[loc_i]->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locClonedFinalParticles[loc_i]);
		else if(locClonedFinalParticles[loc_i]->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locClonedFinalParticles[loc_i], true));
			else
				locNoConstrainParticles.push_back(locClonedFinalParticles[loc_i]);
		}
		else if(locClonedFinalParticles[loc_i]->Get_Charge() == 0)
			locOnlyConstrainTimeParticles.push_back(locClonedFinalParticles[loc_i]);
		else
			locConstrainSpacetimeParticles.push_back(locClonedFinalParticles[loc_i]);
	}

	//set momentum of neutral showers that have 5x5 covariance matrix (shower energy input instead of p3)
	for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	{
		if(!locNoConstrainParticles[loc_i]->Get_IsNeutralShowerFlag())
			continue; //only do for neutral showers
		double locE = locNoConstrainParticles[loc_i]->Get_ShowerEnergy();
		double locMass = locNoConstrainParticles[loc_i]->Get_Mass();
		double locPMag = sqrt(locE*locE - locMass*locMass);
		TVector3 locMomentum = locNoConstrainParticles[loc_i]->Get_Position() - locNoConstrainParticles[loc_i]->Get_CommonVertex();
		locMomentum.SetMag(locPMag);
		locNoConstrainParticles[loc_i]->Set_Momentum(locMomentum);
	}

	//create the constraint and set its members
	DKinFitConstraint_Spacetime* locKinFitConstraint = Get_KinFitConstraintSpacetimeResource();
	locKinFitConstraint->Set_ConstrainSpacetimeParticles(locConstrainSpacetimeParticles);
	locKinFitConstraint->Set_OnlyConstrainTimeParticles(locOnlyConstrainTimeParticles);
	locKinFitConstraint->Set_DecayingParticles(locDecayingParticles);
	locKinFitConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locKinFitConstraint->Set_CommonVertex(locVertexGuess);
	locKinFitConstraint->Set_CommonTime(locCommonTimeGuess);
	locKinFitConstraint->Set_UseRFTimeFlag(locUseRFTimeFlag);

	//add constraint to particles
	TVector3 locMomentum;
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
		locClonedInitialParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locKinFitConstraint);
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
		locClonedFinalParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locKinFitConstraint);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Spacetime constraint added. Vertex/Time constrained particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locConstrainSpacetimeParticles.size(); ++loc_i)
	 		cout << locConstrainSpacetimeParticles[loc_i]->Get_Charge() << ", " << locConstrainSpacetimeParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Time-only constrained particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locOnlyConstrainTimeParticles.size(); ++loc_i)
	 		cout << locOnlyConstrainTimeParticles[loc_i]->Get_Charge() << ", " << locOnlyConstrainTimeParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	 		cout << locNoConstrainParticles[loc_i]->Get_Charge() << ", " << locNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	 		cout << locDecayingParticles[loc_i].first->Get_Charge() << ", " << locDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locKinFitConstraint));
	return locKinFitConstraint;
}

const DKinFitConstraint_P4* DKinFitter::Add_P4Constraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles)
{
	//require p4 is conserved between the tracks

	//first check to make sure the inputs are ok //can't tell if constrained properly until Resolve_Constraints() //due to decaying particles in 2 constraints
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		if(locInitialParticles[loc_i] == NULL)
			return NULL;
	}
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		if(locFinalParticles[loc_i] == NULL)
			return NULL;
	}

	//constraint will be created.  first clone the particles (if not previously cloned)
	deque<DKinFitParticle*> locClonedInitialParticles, locClonedFinalParticles;
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
		locClonedInitialParticles.push_back(GetOrCreate_ClonedParticle(locInitialParticles[loc_i]));
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
		locClonedFinalParticles.push_back(GetOrCreate_ClonedParticle(locFinalParticles[loc_i]));

	//create the constraint and set its members
	DKinFitConstraint_P4* locKinFitConstraint = Get_KinFitConstraintP4Resource();
	locKinFitConstraint->Set_InitialParticles(locClonedInitialParticles);
	locKinFitConstraint->Set_FinalParticles(locClonedFinalParticles);

	//mark constraint in particles
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
		locClonedInitialParticles[loc_i]->Set_IsInP4FitFlag(true);
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
		locClonedFinalParticles[loc_i]->Set_IsInP4FitFlag(true);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: P4 constraint added. Initial-state particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	 		cout << locClonedInitialParticles[loc_i]->Get_Charge() << ", " << locClonedInitialParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Final-state particle q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	 		cout << locClonedFinalParticles[loc_i]->Get_Charge() << ", " << locClonedFinalParticles[loc_i]->Get_Mass() << endl;
	}

	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locKinFitConstraint));
	return locKinFitConstraint;
}

void DKinFitter::Set_Particle(DKinFitParticle* locKinFitParticle)
{
	for(size_t loc_j = 0; loc_j < dKinFitParticles.size(); ++loc_j)
	{
		if(locKinFitParticle == dKinFitParticles[loc_j])
			return; //already added
	}
	dKinFitParticles.push_back(locKinFitParticle);
}

bool DKinFitter::Fit_Reaction(void)
{
	if(!Resolve_Constraints())
		return false;

	Set_MatrixSizes();
	Setup_Matrices();
	Fill_InputMatrices();

	double locPreviousChiSq;

	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dEta: " << endl;
		Print_Matrix(dEta);
		cout << "DKinFitter: dVY: " << endl;
		Print_Matrix(dVY);
		cout << "DKinFitter: dXi: " << endl;
		Print_Matrix(dXi);
		for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
			Print_ParticleParams(dKinFitParticles[loc_i]);
	}

	dChiSq = 9.9E99;
	int locNumIterations = -1;
	do
	{
		++locNumIterations;
		if(locNumIterations >= int(dMaxNumIterations))
		{
			if(dDebugLevel > 10)
				cout << "DKinFitter: Exceeded maximum number of iterations.  Returning false." << endl;
			return false; //diverging!
		}

		locPreviousChiSq = dChiSq;
		if(dDebugLevel > 20)
			cout << "DKinFitter: Begin iteration" << endl;

		Calc_dF();
		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dF: " << endl;
			Print_Matrix(dF);
			cout << "DKinFitter: dFdXi: " << endl;
			Print_Matrix(dFdXi);
			cout << "DKinFitter: dFdEta: " << endl;
			Print_Matrix(dFdEta);
		}

		dR = dF + dFdEta*(dY - dEta);

		if(!Calc_dS())
		{
			if(dDebugLevel > 10)
				cout << "DKinFitter: Failed S-matrix inversion. Returning false." << endl;
			return false; // matrix is not invertible
		}

		if(dNumXi > 0)
		{
			if(!Calc_dVXi())
			{
				if(dDebugLevel > 10)
					cout << "DKinFitter: Failed VXi-matrix inversion. Returning false." << endl;
				return false; // matrix is not invertible
			}

			TMatrixD locDeltaXi = -1.0*(*dVXi)*dFdXi_T*dS_Inverse*dR;
			if(dDebugLevel > 20)
			{
				cout << "DKinFitter: locDeltaXi: " << endl;
				Print_Matrix(locDeltaXi);
			}

			dXi += locDeltaXi;
			if(dDebugLevel > 20)
			{
				cout << "DKinFitter: dXi: " << endl;
				Print_Matrix(dXi);
			}

			dLambda = dS_Inverse*(dR + dFdXi*locDeltaXi);
		}
		else
		{
			dLambda = dS_Inverse*dR;
		}

		dLambda_T.Transpose(dLambda);
		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dLambda: " << endl;
			Print_Matrix(dLambda);
		}

		dEta = dY - dVY*dFdEta_T*dLambda;
		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dEta: " << endl;
			Print_Matrix(dEta);
		}

		TMatrixDSym locTempMatrix5 = dS;
		dChiSq = (locTempMatrix5.SimilarityT(dLambda) + 2.0*dLambda_T*dF)(0, 0); //exercise 10.17 (# unknowns < # measureables, so this is faster!)
		if(dDebugLevel > 20)
			cout << "DKinFitter: dChiSq = " << dChiSq << endl;

		Update_ParticleParams(); //input eta & xi info into particle objects
		if(dDebugLevel > 20)
		{
			for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
				Print_ParticleParams(dKinFitParticles[loc_i]);
		}
	}
	while(fabs(dChiSq - locPreviousChiSq) > 0.001);

	if(dChiSq < 0.0)
	{
		if(dDebugLevel > 10)
			cout << "DKinFitter: Negative dChiSq (" << dChiSq << "). Returning false." << endl;
		return false;
	}

	//calculate final uncertainty matrix
	if(dNumXi > 0)
	{
		TMatrixD locH = dFdEta_T*dS_Inverse*dFdXi;
		TMatrixDSym locTempMatrix10 = dS_Inverse;
		TMatrixDSym locTempMatrix11 = *dVXi;
		*dVEta = dVY - (locTempMatrix10.SimilarityT(dFdEta) - locTempMatrix11.Similarity(locH)).Similarity(dVY);
	}
	else
	{
		TMatrixDSym locTempMatrix10 = dS_Inverse;
		*dVEta = dVY - (locTempMatrix10.SimilarityT(dFdEta)).Similarity(dVY);
	}

	dEpsilon = dY - dEta;

	Calc_Pulls();
	dNDF = dNumF - dNumXi;
	dConfidenceLevel = TMath::Prob(dChiSq, dNDF);

	Set_FinalTrackInfo();

	if(dDebugLevel > 5)
		cout << "DKinFitter: Final dChiSq, dNDF, dConfidenceLevel = " << dChiSq << ", " << dNDF << ", " << dConfidenceLevel << endl;

	return true;
}

bool DKinFitter::Calc_dS(void)
{
	TMatrixDSym locTempMatrix = dVY;
	locTempMatrix.Similarity(dFdEta);
	dS = locTempMatrix;
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dS: " << endl;
		Print_Matrix(dS);
		cout << "determinant magnitude = " << fabs(dS.Determinant()) << endl;
	}
	TDecompLU locDecompLU_S(dS);
	//check to make sure that the matrix is decomposable and has a non-zero determinant
	if((!locDecompLU_S.Decompose()) || (fabs(dS.Determinant()) < 1.0E-300))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitter: dS not invertible.  Returning false." << endl;
		return false; // matrix is not invertible
	}
	dS_Inverse = dS;
	dS_Inverse.Invert();
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dS_Inverse: " << endl;
		Print_Matrix(dS_Inverse);
	}

	return true;
}

bool DKinFitter::Calc_dVXi(void)
{
	TMatrixDSym locTempMatrix = dS_Inverse;
	locTempMatrix.SimilarityT(dFdXi);
	dVXi_Inverse = locTempMatrix;
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dVXi_Inverse: " << endl;
		Print_Matrix(dVXi_Inverse);
		cout << "determinant magnitude = " << fabs(dVXi_Inverse.Determinant()) << endl;
	}
	TDecompLU locDecompLU_VXiInv(dVXi_Inverse);
	//check to make sure that the matrix is decomposable and has a non-zero determinant
	if((!locDecompLU_VXiInv.Decompose()) || (fabs(dVXi_Inverse.Determinant()) < 1.0E-300))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitter: dVXi_Inverse not invertible.  Returning false." << endl;
		return false; // matrix is not invertible
	}
	*dVXi = dVXi_Inverse;
	dVXi->Invert();
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dVXi: " << endl;
		Print_Matrix(*dVXi);
	}
	return true;
}

void DKinFitter::Print_Matrix(const TMatrixD& locMatrix) const
{
	for(int loc_i = 0; loc_i < locMatrix.GetNrows(); ++loc_i)
	{
		for(int loc_j = 0; loc_j < locMatrix.GetNcols(); ++loc_j)
			cout << locMatrix(loc_i, loc_j) << ", ";
		cout << endl;
	}
}

void DKinFitter::Print_ParticleParams(const DKinFitParticle* locKinFitParticle) const
{
	int locCharge = locKinFitParticle->Get_Charge();
	double locMass = locKinFitParticle->Get_Mass();
	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TLorentzVector locSpacetimeVertex = locKinFitParticle->Get_SpacetimeVertex();
	const TMatrixDSym* locCovarianceMatrix = locKinFitParticle->Get_CovarianceMatrix();
	cout << "DKinFitter: Particle Q, Mass, E, P3, V3, T = " << locCharge << ", " << locMass << ", " << locP4.E() << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;
	if(locCovarianceMatrix != NULL)
	{
		cout << "DKinFitter: CovMatrix Diagonal Terms: ";
		for(int loc_i = 0; loc_i < locCovarianceMatrix->GetNcols(); ++loc_i)
			cout << (*locCovarianceMatrix)(loc_i, loc_i) << ", ";
		cout << endl;
	}
	cout << "DKinFitter: Particle E, Px, Vx, T, L indices = " << locKinFitParticle->Get_EParamIndex() << ", " << locKinFitParticle->Get_PxParamIndex() << ", " << locKinFitParticle->Get_VxParamIndex() << ", " << locKinFitParticle->Get_TParamIndex() << ", " << locKinFitParticle->Get_LParamIndex() << endl;
	cout << "DKinFitter: Particle CovMatrix E, Px, Vx, T indices = " << locKinFitParticle->Get_CovMatrixEParamIndex() << ", " << locKinFitParticle->Get_CovMatrixPxParamIndex() << ", " << locKinFitParticle->Get_CovMatrixVxParamIndex() << ", " << locKinFitParticle->Get_CovMatrixTParamIndex() << endl;
}

bool DKinFitter::Resolve_Constraints(void)
{
	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	//make sure there are no neutral showers used in p4 constraints that are not included in a vertex fit
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		if(!locKinFitParticle->Get_IsNeutralShowerFlag())
			continue; //is non-shower
		if(!locKinFitParticle->Get_IsInP4FitFlag())
			continue; //is not in p4 fit
		if(!locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag())
		{
			cout << "ERROR in DKinFitter: Detected neutral shower in a P4 fit but not in a vertex or spacetime fit: P3 is undefined!!  Exiting." << endl;
			return false; //detected neutral shower in a p4 fit but not in a vertex fit: p3 is undefined!!
		}
	}

	//some checks on whether particles setup correctly (not too many missing particles), decaying particle used in vertex constraints correctly
	deque<DKinFitParticle*> locUnconstrainedParticles;
	deque<DKinFitParticle*> locMissingParticles;
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
			continue;

		if(!locKinFitParticle->Get_IsInP4FitFlag())
		{
			//not constrained by p4 fit
			if(locKinFitParticleType == d_MissingParticle)
				continue; //if it's missing but not in a p4 constraint then it's ok: user could just be setting it's vertex or time to the kinfit result (e.g. resonance)
			if((locKinFitParticle->Get_CommonVertexAndOrTimeConstraints().size() < 2) || (!dLinkVerticesFlag))
				continue; //if it's decaying but not in a p4 constraint then it's ok as long as not in 2+ vertex constraints: user could just be setting it's vertex or time to the kinfit result (e.g. resonance)
			cout << "ERROR in DKinFitter: Decaying particle constrained in a vertex or vertex-time fit with unconstrained momentum!!  Exiting." << endl;
			return false; //decaying but constrained in a vertex or vertex-time fit with unconstrained momentum
		}
		locUnconstrainedParticles.push_back(locKinFitParticle);
		if(locKinFitParticleType == d_MissingParticle)
			locMissingParticles.push_back(locKinFitParticle);
	}
	if(locMissingParticles.size() > 1)
	{
		cout << "ERROR in DKinFitter: Cannot constrain more than one missing particle!!  Exiting." << endl;
		return false; //can't constrain more than one missing particle!
	}

	//loop through the p4 constraints, find constraints which contain only one missing or decaying particle: used as starting point for assigning them to different constraints
	deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> > locConstrainableParticles;
	deque<const DKinFitParticle*> locConstrainedParticles;
	while(locConstrainedParticles.size() < locUnconstrainedParticles.size())
	{
		if(locConstrainableParticles.empty())
		{
			if(!Find_ConstrainableParticles(locConstrainableParticles, locConstrainedParticles))
				return false; //more decaying/missing particles than constraints, or circular dependency: cannot fit
		}

		while(!locConstrainableParticles.empty())
		{
			bool locAlreadyConstrainedFlag = false;
			//see if the particle is already constrained (e.g. constrainable in more than one constraint)
			for(size_t loc_j = 0; loc_j < locConstrainedParticles.size(); ++loc_j)
			{
				if(locConstrainableParticles.back().first != locConstrainedParticles[loc_j])
					continue;
				locAlreadyConstrainedFlag = true;
				break;
			}
			if(locAlreadyConstrainedFlag)
			{
				locConstrainableParticles.pop_back();
				continue;
			}

			Constrain_Particle(locConstrainableParticles.back());
			locConstrainedParticles.push_back(locConstrainableParticles.back().first);
			locConstrainableParticles.pop_back();
		}
	}

	if(!Resolve_DecayingParticleSpacetimeLinks())
		return false;

	return true;
}

bool DKinFitter::Find_ConstrainableParticles(deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles)
{
	DKinFitParticleType locKinFitParticleType;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 == NULL)
			continue;
		if(locKinFitConstraint_P4->Get_ConstrainedP4Particle() != NULL)
			continue; //constraint already set
		size_t locNumParticlesNeedToBeConstrained = 0;
		DKinFitParticle* locConstrainableParticle = NULL;
		deque<DKinFitParticle*> locTempParticles = locKinFitConstraint_P4->Get_InitialParticles();
		for(size_t loc_j = 0; loc_j < locTempParticles.size(); ++loc_j)
		{
			locKinFitParticleType = locTempParticles[loc_j]->Get_KinFitParticleType();
			if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
				continue;
			bool locParticleAlreadyConstrainedFlag = false;
			for(size_t loc_k = 0; loc_k < locConstrainedParticles.size(); ++loc_k)
			{
				if(locConstrainedParticles[loc_k] != locTempParticles[loc_j])
					continue;
				locParticleAlreadyConstrainedFlag = true;
				break;
			}
			if(locParticleAlreadyConstrainedFlag)
				continue;
			++locNumParticlesNeedToBeConstrained;
			locConstrainableParticle = locTempParticles[loc_j];
		}
		locTempParticles = locKinFitConstraint_P4->Get_FinalParticles();
		for(size_t loc_j = 0; loc_j < locTempParticles.size(); ++loc_j)
		{
			locKinFitParticleType = locTempParticles[loc_j]->Get_KinFitParticleType();
			if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
				continue;
			bool locParticleAlreadyConstrainedFlag = false;
			for(size_t loc_k = 0; loc_k < locConstrainedParticles.size(); ++loc_k)
			{
				if(locConstrainedParticles[loc_k] != locTempParticles[loc_j])
					continue;
				locParticleAlreadyConstrainedFlag = true;
				break;
			}
			if(locParticleAlreadyConstrainedFlag)
				continue;
			++locNumParticlesNeedToBeConstrained;
			locConstrainableParticle = locTempParticles[loc_j];
		}
		if(locNumParticlesNeedToBeConstrained == 1) //else too many unconstrained particles in it's step to be able to constrain it right away!!
			locConstrainableParticles.push_back(pair<DKinFitParticle*, DKinFitConstraint_P4*>(locConstrainableParticle, locKinFitConstraint_P4));
	}
	return (!locConstrainableParticles.empty());
}

void DKinFitter::Constrain_Particle(pair<DKinFitParticle*, DKinFitConstraint_P4*>& locConstrainableParticle)
{
	DKinFitParticle* locParticleToConstrain = locConstrainableParticle.first;
	DKinFitConstraint_P4* locKinFitConstraint_P4 = locConstrainableParticle.second;
	DKinFitParticle* locKinFitParticle;

	if(dDebugLevel > 5)
		cout << "particle to constrain q, mass = " << locParticleToConstrain->Get_Charge() << ", " << locParticleToConstrain->Get_Mass() << endl;

	TVector3 locMomentum;
	bool locConstrainedParticleIsInInitialState = false;
	for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->Get_InitialParticles().size(); ++loc_j)
	{
		locKinFitParticle = locKinFitConstraint_P4->Get_InitialParticles()[loc_j];
		if(locKinFitParticle == locParticleToConstrain)
		{
			locKinFitConstraint_P4->Set_ConstrainedP4Particle(locParticleToConstrain);
			locParticleToConstrain->Set_DecayingParticleAtProductionVertexFlag(false);
			locConstrainedParticleIsInInitialState = true;
		}
		else
			locMomentum += locKinFitParticle->Get_Momentum();
	}
	for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->Get_FinalParticles().size(); ++loc_j)
	{
		locKinFitParticle = locKinFitConstraint_P4->Get_FinalParticles()[loc_j];
		if(locKinFitParticle == locParticleToConstrain)
		{
			locKinFitConstraint_P4->Set_ConstrainedP4Particle(locParticleToConstrain);
			locParticleToConstrain->Set_DecayingParticleAtProductionVertexFlag(true);
		}
		else
			locMomentum -= locKinFitParticle->Get_Momentum();
	}

	//set p3 guess
	if(locConstrainedParticleIsInInitialState)
		locMomentum *= -1.0;
	locParticleToConstrain->Set_Momentum(locMomentum);
}

bool DKinFitter::Resolve_DecayingParticleSpacetimeLinks(void)
{
	//resolve links between vertex & time fits (decaying particles)
	deque<DKinFitConstraint_Vertex*> locVertexConstraints;
	deque<DKinFitConstraint_Spacetime*> locSpacetimeConstraints;

	//loop over all particles, mostly ignore all but decaying particles
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		//get and sort it's vertex constraints
		deque<DKinFitConstraint_VertexBase*> locVertexAndOrTimeConstraints = dKinFitParticles[loc_i]->Get_CommonVertexAndOrTimeConstraints();
		if(dKinFitParticles[loc_i]->Get_KinFitParticleType() != d_DecayingParticle)
		{
			if(locVertexAndOrTimeConstraints.size() > 1)
			{
				cout << "ERROR in DKinFitter: Detected particle involved in more than 1 vertex constraint!!  Exiting." << endl;
				return false; //detected particle involved in more than 1 vertex constraint
			}
			continue;
		}

		if(locVertexAndOrTimeConstraints.empty())
			continue;
		if(locVertexAndOrTimeConstraints.size() > 2)
		{
			cout << "ERROR in DKinFitter: Decaying particle involved in more than 2 vertex constraints!!  Exiting." << endl;
			return false; //decaying particle involved in more than 2 vertex constraints
		}

		//if the decaying particle is only in one vertex constraint, then it's position is defined by the fit result; it is not used to constrain the fit
		if(locVertexAndOrTimeConstraints.size() == 1)
		{
			DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(locVertexAndOrTimeConstraints[0]);
			if(locKinFitConstraint_Vertex != NULL)
				locKinFitConstraint_Vertex->Add_NoConstrainParticle(dKinFitParticles[loc_i]);
			else //vertex & time fit
				(dynamic_cast<DKinFitConstraint_Spacetime*>(locVertexAndOrTimeConstraints[0]))->Add_NoConstrainParticle(dKinFitParticles[loc_i]);
			continue;
		}

		//make sure it's p4 is constrained since it's in 2 vertex fits
		bool locP4ConstrainedFlag = false;
		for(size_t loc_j = 0; loc_j < dKinFitConstraints.size(); ++loc_j)
		{
			DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_j]);
			if(locKinFitConstraint_P4 == NULL)
				continue;
			if(locKinFitConstraint_P4->Get_ConstrainedP4Particle() != dKinFitParticles[loc_i])
				continue; //either no particle p4 constrained, or is the wrong particle
			locP4ConstrainedFlag = true;
			break;
		}
		if(!locP4ConstrainedFlag)
		{
			cout << "ERROR in DKinFitter: Unable to include decaying particle in 2 vertex/time constraints if it's momentum is unknown (not constrained)!!  Exiting." << endl;
			return false; //unable to include decaying particle in vertex/time constraint if it's momentum is unknown (not constrained)
		}

		//below flag true if decaying particle is an initial particle in the p4 constraint it's set in (defined at)
		bool locInitialStateFlag_P4 = !(dKinFitParticles[loc_i]->Get_DecayingParticleAtProductionVertexFlag());

		//loop over the vertex constraints that the decaying particle is involved in
		bool locKnownVertexNotDefinedFlag = true;
		deque<DKinFitConstraint_VertexBase*>::iterator locCommonVertexConstraintIterator = locVertexAndOrTimeConstraints.end();
		for(deque<DKinFitConstraint_VertexBase*>::iterator locIterator = locVertexAndOrTimeConstraints.begin(); locIterator != locVertexAndOrTimeConstraints.end(); ++locIterator)
		{
			DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(*locIterator);
			if(locKinFitConstraint_Vertex != NULL)
			{
				//below flag true if decaying particle is an initial particle in the vertex constraint
				bool locInitialStateFlag_Vertex = locKinFitConstraint_Vertex->Get_DecayingParticleInInitialStateFlag(dKinFitParticles[loc_i]);

				//the particle position/time is DEFINED at the location the momentum is defined, so the common vertex/time constraint should be on the OTHER position/time of the track
				if(locInitialStateFlag_P4 == locInitialStateFlag_Vertex)
				{
					//the decaying particle is in the same side of the reaction in both constraints
					locKinFitConstraint_Vertex->Add_NoConstrainParticle(dKinFitParticles[loc_i]);
					if(dDebugLevel > 10)
						cout << "q, mass of decaying particle vertex defined by constraint: " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
					locKnownVertexNotDefinedFlag = false;
				}
				else
				{
					locKinFitConstraint_Vertex->Add_ConstrainVertexParticle(dKinFitParticles[loc_i]);
					if(dDebugLevel > 10)
						cout << "q, mass of decaying particle constrained to vertex: " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
					locCommonVertexConstraintIterator = locIterator;
				}
				continue;
			}
			DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(*locIterator);
			if(locKinFitConstraint_Spacetime != NULL)
			{
				//below flag true if decaying particle is an initial particle in the vertex constraint
				bool locInitialStateFlag_Vertex = locKinFitConstraint_Spacetime->Get_DecayingParticleInInitialStateFlag(dKinFitParticles[loc_i]);

				//the particle position/time is DEFINED at the location the momentum is defined, so the common vertex/time constraint should be on the OTHER position/time of the track
				if(locInitialStateFlag_P4 == locInitialStateFlag_Vertex)
				{
					//the decaying particle is in the same side of the reaction in both constraints
					locKinFitConstraint_Spacetime->Add_NoConstrainParticle(dKinFitParticles[loc_i]);
					if(dDebugLevel > 10)
						cout << "q, mass of decaying particle spacetime vertex defined by constraint: " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
					locKnownVertexNotDefinedFlag = false;
				}
				else
				{
					locKinFitConstraint_Spacetime->Add_ConstrainSpacetimeParticle(dKinFitParticles[loc_i]);
					if(dDebugLevel > 10)
						cout << "q, mass of decaying particle constrained to spacetime vertex: " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
					locCommonVertexConstraintIterator = locIterator;
				}
				continue;
			}
		}
		if(locKnownVertexNotDefinedFlag)
		{
			cout << "ERROR in DKinFitter: decaying particle constrained to a vertex, but no point on its trajectory is defined (by another vertex constraint)!!  Exiting." << endl;
			return false; //decaying particle constrained to a vertex, but no point on its trajectory is defined (by another vertex constraint)
		}

		//sort the vertex constraints
		DKinFitConstraint_VertexBase* locTempConstraint = *locCommonVertexConstraintIterator;
		locVertexAndOrTimeConstraints.erase(locCommonVertexConstraintIterator);
		locVertexAndOrTimeConstraints.push_front(locTempConstraint);
		dKinFitParticles[loc_i]->Set_CommonVertexAndOrTimeConstraints(locVertexAndOrTimeConstraints);
	}

	//make sure that each vertex has at least 2 particles to constrain it
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			if(locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles().size() < 2)
			{
				cout << "ERROR in DKinFitter: not enough particles in spacetime constraint!!  Exiting." << endl;
				return false; //decaying particle constrained to a vertex, but no point on its trajectory is defined (by another vertex constraint)
			}
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			if(locKinFitConstraint_Vertex->Get_ConstrainVertexParticles().size() < 2)
			{
				cout << "ERROR in DKinFitter: not enough particles in spacetime constraint!!  Exiting." << endl;
				return false; //decaying particle constrained to a vertex, but no point on its trajectory is defined (by another vertex constraint)
			}
		}
	}

	return true;
}

void DKinFitter::Set_MatrixSizes(void)
{
	//set matrix sizes
	dNumXi = 0; //num unknowns
	dNumEta = 0; //num measurables
	dNumF = 0; //num constraint eqs

	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	//Calculate dNumEta
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		if(!locKinFitParticle->Get_IsNeutralShowerFlag())
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) || (locKinFitParticle->Get_IsConstrainedByVertexOrSpacetimeFitFlag()))
				dNumEta += 3; //p3
			if(locKinFitParticle->Get_IsConstrainedByVertexOrSpacetimeFitFlag())
				dNumEta += 3; //v3
		}
		else //neutral shower
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) && (locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag()))
				dNumEta += 4; //E + v3 (p4 fit needs p3, which is derived from v3 + kinfit vertex)
		}
		if(locKinFitParticle->Get_IsInSpacetimeFitFlag())
			++dNumEta; //t
	}
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i) //check if RF time
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime == NULL)
			continue;
		if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
			++dNumEta; //RF t
	}

	//Calculate dNumXi and dNumF
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			dNumF += 4; //p4
			if(locKinFitConstraint_P4->Get_ConstrainedP4Particle() != NULL)
				dNumXi += 3; //p3 //missing particle
			continue;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			dNumXi += 3; //v3
			dNumF += 2*(locKinFitConstraint_Vertex->Get_ConstrainVertexParticles().size()); //for each track
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->Get_ConstrainVertexParticles().size(); ++loc_j)
					cout << locKinFitConstraint_Vertex->Get_ConstrainVertexParticles()[loc_j]->Get_Charge() << ", " << locKinFitConstraint_Vertex->Get_ConstrainVertexParticles()[loc_j]->Get_Mass() << ";";
				cout << endl;
			}
			continue;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			dNumXi += 4; //v3, t
			dNumF += 3*(locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles().size()); //for each track
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of spacetime vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles().size(); ++loc_j)
					cout << locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles()[loc_j]->Get_Charge() << ", " << locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles()[loc_j]->Get_Mass() << ";";
				cout << endl;
			}
			if(Get_IsBFieldNearBeamline())
			{
				size_t locNumChargedConstraintParticles = 0;
				deque<DKinFitParticle*> locConstrainSpacetimeParticles = locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles();
				for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
				{
					if(locConstrainSpacetimeParticles[loc_j]->Get_Charge() != 0)
						++locNumChargedConstraintParticles;
				}
				dNumXi += locNumChargedConstraintParticles; //path length (l) for each accelerating particle
				dNumF += locNumChargedConstraintParticles; //extra constraint due to extra unknown (path length (l) for each accelerating particle)
			}
			dNumF += locKinFitConstraint_Spacetime->Get_OnlyConstrainTimeParticles().size(); //for each neutral shower
			if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
				++dNumF;
		}
	}

	if(dDebugLevel > 10)
		cout << "DKinFitter: Num measurables, unknowns, constraints = " << dNumEta << ", " << dNumXi << ", " << dNumF << endl;
}

void DKinFitter::Setup_Matrices(void)
{
	if(dF.GetNrows() != static_cast<int>(dNumF))
	{
		dF.ResizeTo(dNumF, 1);
		dR.ResizeTo(dNumF, 1);
		dS.ResizeTo(dNumF, dNumF);
		dS_Inverse.ResizeTo(dNumF, dNumF);
		dLambda.ResizeTo(dNumF, 1);
		dLambda_T.ResizeTo(1, dNumF);
		dFdEta.ResizeTo(dNumF, dNumEta);
		dFdEta_T.ResizeTo(dNumEta, dNumF);
		dFdXi.ResizeTo(dNumF, dNumXi);
		dFdXi_T.ResizeTo(dNumXi, dNumF);
	}
	else
	{
		if(dFdEta.GetNcols() != static_cast<int>(dNumEta))
		{
			dFdEta.ResizeTo(dNumF, dNumEta);
			dFdEta_T.ResizeTo(dNumEta, dNumF);
		}
		if(dFdXi.GetNcols() != static_cast<int>(dNumXi))
		{
			dFdXi.ResizeTo(dNumF, dNumXi);
			dFdXi_T.ResizeTo(dNumXi, dNumF);
		}
	}

	if(dY.GetNrows() != static_cast<int>(dNumEta))
	{
		dY.ResizeTo(dNumEta, 1);
		dEta.ResizeTo(dNumEta, 1);
		dEpsilon.ResizeTo(dNumEta, 1);
		dVY.ResizeTo(dNumEta, dNumEta);
	}
	dVEta->ResizeTo(dNumEta, dNumEta);

	if(dXi.GetNrows() != static_cast<int>(dNumXi))
	{
		dXi.ResizeTo(dNumXi, 1);
		dVXi_Inverse.ResizeTo(dNumXi, dNumXi);
	}
	dVXi->ResizeTo(dNumXi, dNumXi);

	Zero_Matrices(); //zeroes all class matrices
}

void DKinFitter::Zero_Matrices(void)
{
	dXi.Zero();
	dEta.Zero();
	dY.Zero();
	dVY.Zero();
	dF.Zero();
	dEpsilon.Zero();

	dLambda.Zero();
	dLambda_T.Zero();
	dFdEta.Zero();
	dFdEta_T.Zero();
	dFdXi.Zero();
	dFdXi_T.Zero();

	dS.Zero();
	dS_Inverse.Zero();
	dR.Zero();

	dVXi->Zero();
	dVXi_Inverse.Zero();
	dVEta->Zero();
}

void DKinFitter::Fill_InputMatrices(void)
{
	//fill dY, dEta, dVY, dXi

	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	int locCharge, locParamIndex, locConstraintIndex;
	TVector3 locMomentum, locPosition;

	//SETUP dY
	locParamIndex = 0;
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		locMomentum = locKinFitParticle->Get_Momentum();
		locPosition = locKinFitParticle->Get_Position();
		locCharge = locKinFitParticle->Get_Charge();

		if(!locKinFitParticle->Get_IsNeutralShowerFlag()) //non-neutral shower
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) || (locKinFitParticle->Get_IsConstrainedByVertexOrSpacetimeFitFlag())) //p3
			{
				locKinFitParticle->Set_PxParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locMomentum.Px();
				dY(locParamIndex + 1, 0) = locMomentum.Py();
				dY(locParamIndex + 2, 0) = locMomentum.Pz();
				locParamIndex += 3;
			}
			if(locKinFitParticle->Get_IsConstrainedByVertexOrSpacetimeFitFlag())
			{
				locKinFitParticle->Set_VxParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locPosition.Px();
				dY(locParamIndex + 1, 0) = locPosition.Py();
				dY(locParamIndex + 2, 0) = locPosition.Pz();
				locParamIndex += 3;
			}
		}
		else //neutral shower
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) && (locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag()))
			{
				//E + v3 (p4 fit needs p3, which is derived from v3 + kinfit vertex)
				locKinFitParticle->Set_EParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locKinFitParticle->Get_ShowerEnergy();
				++locParamIndex;

				locKinFitParticle->Set_VxParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locPosition.Px();
				dY(locParamIndex + 1, 0) = locPosition.Py();
				dY(locParamIndex + 2, 0) = locPosition.Pz();
				locParamIndex += 3;
			}
		}

		if(locKinFitParticle->Get_IsInSpacetimeFitFlag())
		{
			locKinFitParticle->Set_TParamIndex(locParamIndex);
			dY(locParamIndex, 0) = locKinFitParticle->Get_Time();
			++locParamIndex;
		}
	}
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i) //set RF time
	{
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime == NULL)
			continue;
		if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
		{
			dRFTimeParamIndex = locParamIndex;
			dY(locParamIndex, 0) = dRFTime;
			++locParamIndex;
		}
	}

	//SETUP dEta
	dEta = dY; //use measurements as first guess

	//SETUP dXi (with initial guesses)
	locParamIndex = 0;
	locConstraintIndex = 0;
	DKinFitParticle* locConstrainedKinFitParticle;
	deque<DKinFitParticle*> locNoConstrainParticles;
	bool locRFTimeConstrainedFlag = false;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			locKinFitConstraint_P4->Set_FIndex(locConstraintIndex);
			locConstraintIndex += 4;
			locConstrainedKinFitParticle = locKinFitConstraint_P4->Get_ConstrainedP4Particle();
			if(locConstrainedKinFitParticle == NULL)
				continue; //e.g. no missing particles
			//set initial p3 guess
			TVector3 locMomentum = locConstrainedKinFitParticle->Get_Momentum();
			dXi(locParamIndex, 0) = locMomentum.Px();
			dXi(locParamIndex + 1, 0) = locMomentum.Py();
			dXi(locParamIndex + 2, 0) = locMomentum.Pz();
			locConstrainedKinFitParticle->Set_PxParamIndex(locParamIndex);
			locParamIndex += 3;
			continue;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			locPosition = locKinFitConstraint_Vertex->Get_CommonVertex();
			dXi(locParamIndex, 0) = locPosition.X();
			dXi(locParamIndex + 1, 0) = locPosition.Y();
			dXi(locParamIndex + 2, 0) = locPosition.Z();
			locKinFitConstraint_Vertex->Set_VxParamIndex(locParamIndex);
			locNoConstrainParticles = locKinFitConstraint_Vertex->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
			{
				if((locNoConstrainParticles[loc_j]->Get_KinFitParticleType() == d_MissingParticle) || (locNoConstrainParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle))
					locNoConstrainParticles[loc_j]->Set_VxParamIndex(locParamIndex); //not included in fit, but particle vertex is defined by the fit result
			}
			locParamIndex += 3;
			locKinFitConstraint_Vertex->Set_FIndex(locConstraintIndex);
			locConstraintIndex += 2*(locKinFitConstraint_Vertex->Get_ConstrainVertexParticles().size());
			continue;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			locPosition = locKinFitConstraint_Spacetime->Get_CommonVertex();
			dXi(locParamIndex, 0) = locPosition.X();
			dXi(locParamIndex + 1, 0) = locPosition.Y();
			dXi(locParamIndex + 2, 0) = locPosition.Z();
			dXi(locParamIndex + 3, 0) = locKinFitConstraint_Spacetime->Get_CommonTime();
			locKinFitConstraint_Spacetime->Set_VxParamIndex(locParamIndex);
			locKinFitConstraint_Spacetime->Set_TParamIndex(locParamIndex + 3);

			locNoConstrainParticles = locKinFitConstraint_Vertex->Get_NoConstrainParticles();
			for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
			{
				if((locNoConstrainParticles[loc_j]->Get_KinFitParticleType() != d_MissingParticle) && (locNoConstrainParticles[loc_j]->Get_KinFitParticleType() != d_DecayingParticle))
					continue;
				locNoConstrainParticles[loc_j]->Set_VxParamIndex(locParamIndex); //not included in fit, but particle vertex is defined by the fit result
				locNoConstrainParticles[loc_j]->Set_TParamIndex(locParamIndex + 3); //not included in fit, but particle time is defined by the fit result
			}
			locParamIndex += 4;

			locKinFitConstraint_Spacetime->Set_FIndex(locConstraintIndex);
			locConstraintIndex += 3*(locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles().size());
			if(Get_IsBFieldNearBeamline())
			{
				deque<DKinFitParticle*> locConstrainSpacetimeParticles = locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles();
				for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
				{
					if(locConstrainSpacetimeParticles[loc_j]->Get_Charge() == 0)
						continue;
					locConstrainSpacetimeParticles[loc_j]->Set_LParamIndex(locParamIndex);
					dXi(locParamIndex, 0) = locConstrainSpacetimeParticles[loc_j]->Get_PathLength();
					++locParamIndex;
					++locConstraintIndex;
				}
			}
			locConstraintIndex += locKinFitConstraint_Spacetime->Get_OnlyConstrainTimeParticles().size(); //for each neutral shower
			if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
			{
				locRFTimeConstrainedFlag = true;
				++locConstraintIndex;
			}
			continue;
		}
	}

	//SETUP dVY
	int locPxParamIndex;
	int locVxParamIndex;
	int locTParamIndex;
	int locEParamIndex;

	int locCovMatrixEParamIndex;
	int locCovMatrixPxParamIndex;
	int locCovMatrixVxParamIndex;
	int locCovMatrixTParamIndex;
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue; //uncertainties of target particle momentum are zero

		locMomentum = locKinFitParticle->Get_Momentum();
		locPosition = locKinFitParticle->Get_Position();
		const TMatrixDSym& locCovarianceMatrix = *(locKinFitParticle->Get_CovarianceMatrix());

		locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
		locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
		locTParamIndex = locKinFitParticle->Get_TParamIndex();
		locEParamIndex = locKinFitParticle->Get_EParamIndex();

		locCovMatrixEParamIndex = locKinFitParticle->Get_CovMatrixEParamIndex();
		locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
		locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
		locCovMatrixTParamIndex = locKinFitParticle->Get_CovMatrixTParamIndex();

		//localized terms (E, p, v, t)
		if(locEParamIndex >= 0)
			dVY(locEParamIndex, locEParamIndex) = locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixEParamIndex);
		if(locPxParamIndex >= 0)
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
					dVY(locPxParamIndex + loc_j, locPxParamIndex + loc_k) = locCovarianceMatrix(loc_j + locCovMatrixPxParamIndex, loc_k + locCovMatrixPxParamIndex);
			}
		}
		if(locVxParamIndex >= 0)
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
					dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_k) = locCovarianceMatrix(loc_j + locCovMatrixVxParamIndex, loc_k + locCovMatrixVxParamIndex);
			}
		}
		if(locTParamIndex >= 0)
			dVY(locTParamIndex, locTParamIndex) = locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixTParamIndex);

		//cross terms
		if((locEParamIndex >= 0) && (locVxParamIndex >= 0))
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				dVY(locEParamIndex + 0, locVxParamIndex + loc_j) = locCovarianceMatrix(locCovMatrixEParamIndex + 0, locCovMatrixVxParamIndex + loc_j);
				dVY(locVxParamIndex + loc_j, locEParamIndex + 0) = locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixEParamIndex + 0);
			}
		}
		if((locEParamIndex >= 0) && (locTParamIndex >= 0))
		{
			dVY(locEParamIndex, locTParamIndex) = locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex);
			dVY(locTParamIndex, locEParamIndex) = locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex);
		}
		if((locPxParamIndex >= 0) && (locVxParamIndex >= 0))
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				{
					dVY(locPxParamIndex + loc_j, locVxParamIndex + loc_k) = locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixVxParamIndex + loc_k);
					dVY(locVxParamIndex + loc_k, locPxParamIndex + loc_j) = locCovarianceMatrix(locCovMatrixVxParamIndex + loc_k, locCovMatrixPxParamIndex + loc_j);
				}
			}
		}
		if((locPxParamIndex >= 0) && (locTParamIndex >= 0))
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				dVY(locPxParamIndex + loc_j, locTParamIndex + 0) = locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0);
				dVY(locTParamIndex + 0, locPxParamIndex + loc_j) = locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j);
			}
		}
		if((locVxParamIndex >= 0) && (locTParamIndex >= 0))
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				dVY(locVxParamIndex + loc_j, locTParamIndex + 0) = locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0);
				dVY(locTParamIndex + 0, locVxParamIndex + loc_j) = locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j);
			}
		}
	}
	if(locRFTimeConstrainedFlag)
		dVY(dRFTimeParamIndex, dRFTimeParamIndex) = dRFUncertainty*dRFUncertainty;
}

void DKinFitter::Update_ParticleParams(void)
{
	//translate data from Eta to particles
	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	int locParamIndex;

	// update constraint information
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			locParamIndex = locKinFitConstraint_Vertex->Get_VxParamIndex();
			locKinFitConstraint_Vertex->Set_CommonVertex(TVector3(dXi(locParamIndex, 0), dXi(locParamIndex + 1, 0), dXi(locParamIndex + 2, 0)));
			continue;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			locParamIndex = locKinFitConstraint_Spacetime->Get_TParamIndex();
			locKinFitConstraint_Spacetime->Set_CommonTime(dXi(locParamIndex, 0));
			locParamIndex = locKinFitConstraint_Spacetime->Get_VxParamIndex();
			locKinFitConstraint_Spacetime->Set_CommonVertex(TVector3(dXi(locParamIndex, 0), dXi(locParamIndex + 1, 0), dXi(locParamIndex + 2, 0)));

			if(!Get_IsBFieldNearBeamline())
				continue;

			deque<DKinFitParticle*> locConstrainSpacetimeParticles = locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles();
			for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
			{
				if(locConstrainSpacetimeParticles[loc_j]->Get_Charge() == 0)
					continue;
				locParamIndex = locConstrainSpacetimeParticles[loc_j]->Get_LParamIndex();
				if(locParamIndex > 0)
					locConstrainSpacetimeParticles[loc_j]->Set_PathLength(dXi(locParamIndex, 0));
			}
			continue;
		}
	}

	// update particle information
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if(locKinFitParticleType == d_TargetParticle)
			continue; //nothing to update
		else if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle))
		{
			locParamIndex = locKinFitParticle->Get_PxParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Momentum(TVector3(dXi(locParamIndex, 0), dXi(locParamIndex + 1, 0), dXi(locParamIndex + 2, 0)));
			locParamIndex = locKinFitParticle->Get_VxParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Position(TVector3(dXi(locParamIndex, 0), dXi(locParamIndex + 1, 0), dXi(locParamIndex + 2, 0)));
			locParamIndex = locKinFitParticle->Get_TParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Time(dXi(locParamIndex, 0));
		}
		else //initial or detected
		{
			locParamIndex = locKinFitParticle->Get_PxParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Momentum(TVector3(dEta(locParamIndex, 0), dEta(locParamIndex + 1, 0), dEta(locParamIndex + 2, 0)));
			locParamIndex = locKinFitParticle->Get_VxParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Position(TVector3(dEta(locParamIndex, 0), dEta(locParamIndex + 1, 0), dEta(locParamIndex + 2, 0)));
			locParamIndex = locKinFitParticle->Get_TParamIndex();
			if(locParamIndex >= 0)
				locKinFitParticle->Set_Time(dEta(locParamIndex, 0));
			locParamIndex = locKinFitParticle->Get_EParamIndex();
			if(locParamIndex >= 0) //set momentum also //must be after Vx & common vertex are set
			{
				double locE = dEta(locParamIndex, 0);
				locKinFitParticle->Set_ShowerEnergy(locE);
				double locPMag = sqrt(locE*locE - locKinFitParticle->Get_Mass()*locKinFitParticle->Get_Mass());
				TVector3 locMomentum = locKinFitParticle->Get_Position() - locKinFitParticle->Get_CommonVertex();
				locMomentum.SetMag(locPMag);
				locKinFitParticle->Set_Momentum(locMomentum);
			}
		}
	}

	if(dRFTimeParamIndex > 0)
		dRFTime = dEta(dRFTimeParamIndex, 0);
}

void DKinFitter::Calc_dF(void)
{
	dF.Zero();
	dFdXi.Zero();
	dFdEta.Zero();
	size_t locFIndex = 0;
	DKinFitParticle* locKinFitParticle;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		locFIndex = dKinFitConstraints[loc_i]->Get_FIndex();
		if(dDebugLevel > 10)
			cout << "DKinFitter: F index = " << locFIndex << endl;
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->Get_InitialParticles().size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_P4->Get_InitialParticles())[loc_j];
				Calc_dF_P4(locFIndex, locKinFitParticle, true, locKinFitConstraint_P4->Get_ConstrainedP4Particle());
			}
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->Get_FinalParticles().size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_P4->Get_FinalParticles())[loc_j];
				Calc_dF_P4(locFIndex, locKinFitParticle, false, locKinFitConstraint_P4->Get_ConstrainedP4Particle());
			}
			continue;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->Get_ConstrainVertexParticles().size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Vertex->Get_ConstrainVertexParticles())[loc_j];
				Calc_dF_Vertex(locFIndex, locKinFitParticle);
				locFIndex += 2;
			}
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles().size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Spacetime->Get_ConstrainSpacetimeParticles())[loc_j];
				Calc_dF_Time(locFIndex, locKinFitParticle, false); //locFIndex incremented in this func
				if((locKinFitParticle == dRFMatchedBeamParticle) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
					Calc_dF_Time(locFIndex, dRFMatchedBeamParticle, true); //locFIndex incremented in this func
			}
		}
	}
	dFdEta_T.Transpose(dFdEta);
	dFdXi_T.Transpose(dFdXi);
}

void DKinFitter::Calc_dF_P4(size_t locFIndex, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, const DKinFitParticle* locConstrainedParticle)
{
	//E, px, py, pz
	int locCharge = locKinFitParticle->Get_Charge();
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locBField = Get_BField(locPosition);
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;

	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	bool locChargedBFieldVertexFlag = (locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag() && (locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticle != locConstrainedParticle));
	if(locChargedBFieldVertexFlag && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_DecayingParticle)))
	{
		//fitting vertex of charged track in magnetic field that is not the constrained (decaying or missing) particle: px & py change as function of vertex
		locP4.SetPx(locP4.Px() - locA*locH.Z()*locDeltaX.Y() + locA*locH.Y()*locDeltaX.Z());
		locP4.SetPy(locP4.Py() - locA*locH.X()*locDeltaX.Z() + locA*locH.Z()*locDeltaX.X());
		locP4.SetPz(locP4.Pz() - locA*locH.Y()*locDeltaX.X() + locA*locH.X()*locDeltaX.Y());
	}
	//note: neutral shower momentum already calculated & set: done upon constraint input and by Update_ParticleParams
	double locSignMultiplier = locInitialStateFlag ? 1.0 : -1.0;

	dF(locFIndex, 0) += locSignMultiplier*locP4.E();
	dF(locFIndex + 1, 0) += locSignMultiplier*locP4.Px();
	dF(locFIndex + 2, 0) += locSignMultiplier*locP4.Py();
	dF(locFIndex + 3, 0) += locSignMultiplier*locP4.Pz();

	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 1" << endl;
		return; //target params are fixed: no partial derivatives
	}

	if((locKinFitParticle == locConstrainedParticle) || (locKinFitParticleType == d_MissingParticle) || (!locChargedBFieldVertexFlag && (locKinFitParticleType == d_DecayingParticle)))
	{
		//either constrained particle, missing particle, or is a decaying particle constrained by a different constraint but with p3 independent of position (either because neutral, no b-field, or no vertex fit)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 2" << endl;

		dFdXi(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dFdXi(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dFdXi(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dFdXi(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dFdXi(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dFdXi(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;
	}
	else if(locChargedBFieldVertexFlag && (locKinFitParticleType == d_DetectedParticle))
	{
		//detected charged particle in b-field & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 3" << endl;

		dFdEta(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dFdEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dFdEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dFdEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dFdEta(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dFdEta(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;

		dFdEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locA*locH.Z();
		dFdEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locSignMultiplier*locA*locH.Y();

		dFdEta(locFIndex + 2, locVxParamIndex) = -1.0*locSignMultiplier*locA*locH.Z();
		dFdEta(locFIndex + 2, locVxParamIndex + 2) = locSignMultiplier*locA*locH.X();

		dFdEta(locFIndex + 3, locVxParamIndex) = locSignMultiplier*locA*locH.Y();
		dFdEta(locFIndex + 3, locVxParamIndex + 1) = -1.0*locSignMultiplier*locA*locH.X();

		dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dFdEta(locFIndex + 1, locVxParamIndex + 1);
		dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dFdEta(locFIndex + 1, locVxParamIndex + 2);

		dFdXi(locFIndex + 2, locCommonVxParamIndex) -= dFdEta(locFIndex + 2, locVxParamIndex);
		dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dFdEta(locFIndex + 2, locVxParamIndex + 2);

		dFdXi(locFIndex + 3, locCommonVxParamIndex) -= dFdEta(locFIndex + 3, locVxParamIndex);
		dFdXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dFdEta(locFIndex + 3, locVxParamIndex + 1);
	}
	else if(locChargedBFieldVertexFlag && (locKinFitParticleType == d_DecayingParticle))
	{
		//decaying charged particle in b-field, in vertex fit, and p4 constrained by a different constraint
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 4" << endl;

		dFdXi(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dFdXi(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dFdXi(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dFdXi(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dFdXi(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dFdXi(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;

		dFdXi(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locA*locH.Z();
		dFdXi(locFIndex + 1, locVxParamIndex + 2) = -1.0*locSignMultiplier*locA*locH.Y();

		dFdXi(locFIndex + 2, locVxParamIndex) = -1.0*locSignMultiplier*locA*locH.Z();
		dFdXi(locFIndex + 2, locVxParamIndex + 2) = locSignMultiplier*locA*locH.X();

		dFdXi(locFIndex + 3, locVxParamIndex) = locSignMultiplier*locA*locH.Y();
		dFdXi(locFIndex + 3, locVxParamIndex + 1) = -1.0*locSignMultiplier*locA*locH.X();

		dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dFdXi(locFIndex + 1, locVxParamIndex + 1);
		dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dFdXi(locFIndex + 1, locVxParamIndex + 2);

		dFdXi(locFIndex + 2, locCommonVxParamIndex) -= dFdXi(locFIndex + 2, locVxParamIndex);
		dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dFdXi(locFIndex + 2, locVxParamIndex + 2);

		dFdXi(locFIndex + 3, locCommonVxParamIndex) -= dFdXi(locFIndex + 3, locVxParamIndex);
		dFdXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dFdXi(locFIndex + 3, locVxParamIndex + 1);
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 5" << endl;

		dFdEta(locFIndex, locEParamIndex) = locSignMultiplier;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dFdEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Px();
		dFdEta(locFIndex + 2, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Py();
		dFdEta(locFIndex + 3, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Pz();

		TVector3 locDeltaXOverMagDeltaXSq = locDeltaX*(1.0/locDeltaX.Mag2());

		dFdEta(locFIndex + 1, locVxParamIndex) = locSignMultiplier*locP4.Px()*(locDeltaXOverMagDeltaXSq.X() - 1.0/locDeltaX.X());
		dFdEta(locFIndex + 2, locVxParamIndex + 1) = locSignMultiplier*locP4.Py()*(locDeltaXOverMagDeltaXSq.Y() - 1.0/locDeltaX.Y());
		dFdEta(locFIndex + 3, locVxParamIndex + 2) = locSignMultiplier*locP4.Pz()*(locDeltaXOverMagDeltaXSq.Z() - 1.0/locDeltaX.Z());

		dFdEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Y();
		dFdEta(locFIndex + 1, locVxParamIndex + 2) = locSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Z();

		dFdEta(locFIndex + 2, locVxParamIndex) = locSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.X();
		dFdEta(locFIndex + 2, locVxParamIndex + 2) = locSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.Z();

		dFdEta(locFIndex + 3, locVxParamIndex) = locSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.X();
		dFdEta(locFIndex + 3, locVxParamIndex + 1) = locSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.Y();

		dFdXi(locFIndex + 1, locCommonVxParamIndex) -= dFdEta(locFIndex + 1, locVxParamIndex);
		dFdXi(locFIndex + 2, locCommonVxParamIndex + 1) -= dFdEta(locFIndex + 2, locVxParamIndex + 1);
		dFdXi(locFIndex + 3, locCommonVxParamIndex + 2) -= dFdEta(locFIndex + 3, locVxParamIndex + 2);

		dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dFdEta(locFIndex + 1, locVxParamIndex + 1);
		dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dFdEta(locFIndex + 1, locVxParamIndex + 2);

		dFdXi(locFIndex + 2, locCommonVxParamIndex) -= dFdEta(locFIndex + 2, locVxParamIndex);
		dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dFdEta(locFIndex + 2, locVxParamIndex + 2);

		dFdXi(locFIndex + 3, locCommonVxParamIndex) -= dFdEta(locFIndex + 3, locVxParamIndex);
		dFdXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dFdEta(locFIndex + 3, locVxParamIndex + 1);
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 6" << endl;

		dFdEta(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dFdEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dFdEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dFdEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dFdEta(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dFdEta(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;
	}
}

void DKinFitter::Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle)
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
	int locCharge = locKinFitParticle->Get_Charge();

	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //charged track in magnetic field
	{
		TVector3 locBField = Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		TVector3 locPCrossH = locMomentum.Cross(locH);
		TVector3 locPCrossDeltaX = locMomentum.Cross(locDeltaX);
		double locPCrossHMagSq = locPCrossH.Mag2();

		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);
		double locPDotDeltaX = locMomentum.Dot(locDeltaX);

		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
		double locK = locPDotDeltaX - locPDotH*locDeltaXDotH;
		double locJ = locA*locK/locPCrossHMagSq;

		dF(locFIndex, 0) = locPCrossDeltaX.Dot(locH) - 0.5*locA*(locDeltaX.Mag2() - locDeltaXDotH*locDeltaXDotH);
		dF(locFIndex + 1, 0) = locDeltaXDotH - (locPDotH/locA)*asin(locJ);

		TVector3 locM = locDeltaX - locDeltaXDotH*locH;
		double locC = locPDotH/(locPCrossHMagSq*sqrt(1.0 - locJ*locJ));

		TVector3 locdf1dx0 = locPCrossH + locA*locM;
		TVector3 locdf2dx0 = locC*(locMomentum - locPDotH*locH) - locH;

		TVector3 locPCrossHCrossH = locPCrossH.Cross(locH);

		TVector3 locdf2dp0 = -1.0*locH*(asin(locJ)/locA) - locC*(locM + 2.0*(locK/locPCrossHMagSq)*locPCrossHCrossH);

		if(locKinFitParticleType == d_DecayingParticle) //decaying charged particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 1" << endl;

			dFdXi(locFIndex, locPxParamIndex) = locDeltaXCrossH.X();
			dFdXi(locFIndex, locPxParamIndex + 1) = locDeltaXCrossH.Y();
			dFdXi(locFIndex, locPxParamIndex + 2) = locDeltaXCrossH.Z();

			dFdXi(locFIndex, locVxParamIndex) = locdf1dx0.X();
			dFdXi(locFIndex, locVxParamIndex + 1) = locdf1dx0.Y();
			dFdXi(locFIndex, locVxParamIndex + 2) = locdf1dx0.Z();

			dFdXi(locFIndex + 1, locPxParamIndex) = locdf2dp0.X();
			dFdXi(locFIndex + 1, locPxParamIndex + 1) = locdf2dp0.Y();
			dFdXi(locFIndex + 1, locPxParamIndex + 2) = locdf2dp0.Z();

			dFdXi(locFIndex + 1, locVxParamIndex) = locdf2dx0.X();
			dFdXi(locFIndex + 1, locVxParamIndex + 1) = locdf2dx0.Y();
			dFdXi(locFIndex + 1, locVxParamIndex + 2) = locdf2dx0.Z();

			dFdXi(locFIndex, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex, locVxParamIndex);
			dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex, locVxParamIndex + 1);
			dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex, locVxParamIndex + 2);

			dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 1);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 2);
		}
		else //detected particle or beam particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 2" << endl;

			dFdEta(locFIndex, locPxParamIndex) = locDeltaXCrossH.X();
			dFdEta(locFIndex, locPxParamIndex + 1) = locDeltaXCrossH.Y();
			dFdEta(locFIndex, locPxParamIndex + 2) = locDeltaXCrossH.Z();

			dFdEta(locFIndex, locVxParamIndex) = locdf1dx0.X();
			dFdEta(locFIndex, locVxParamIndex + 1) = locdf1dx0.Y();
			dFdEta(locFIndex, locVxParamIndex + 2) = locdf1dx0.Z();

			dFdEta(locFIndex + 1, locPxParamIndex) = locdf2dp0.X();
			dFdEta(locFIndex + 1, locPxParamIndex + 1) = locdf2dp0.Y();
			dFdEta(locFIndex + 1, locPxParamIndex + 2) = locdf2dp0.Z();

			dFdEta(locFIndex + 1, locVxParamIndex) = locdf2dx0.X();
			dFdEta(locFIndex + 1, locVxParamIndex + 1) = locdf2dx0.Y();
			dFdEta(locFIndex + 1, locVxParamIndex + 2) = locdf2dx0.Z();

			dFdXi(locFIndex, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex, locVxParamIndex);
			dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex, locVxParamIndex + 1);
			dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex, locVxParamIndex + 2);

			dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 1);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 2);
		}
	}
	else //neutral particle (e.g. beam photon) or no magnetic field
	{
		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();

			if(locKinFitParticleType == d_DecayingParticle) //decaying particle
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 3" << endl;

				dFdXi(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
				dFdXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
				dFdXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
				dFdXi(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X();

				dFdXi(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
				dFdXi(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
				dFdXi(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
				dFdXi(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();

				dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex, locVxParamIndex + 1);
				dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex, locVxParamIndex + 2);
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 2);
			}
			else
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 4" << endl;

				dFdEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
				dFdEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y(); //soon 0
				dFdEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
				dFdEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X(); //soon 0

				dFdEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
				dFdEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y(); //0
				dFdEta(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
				dFdEta(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X(); //0

				dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex, locVxParamIndex + 1);
				dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex, locVxParamIndex + 2); //0
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 2); //0
			}

		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(locKinFitParticleType == d_DecayingParticle) //decaying particle
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 5" << endl;

				dFdXi(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
				dFdXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
				dFdXi(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
				dFdXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

				dFdXi(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
				dFdXi(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
				dFdXi(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
				dFdXi(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

				dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex, locVxParamIndex + 1);
				dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex, locVxParamIndex + 2);
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 1);
			}
			else
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 6" << endl;

				dFdEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
				dFdEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
				dFdEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
				dFdEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

				dFdEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
				dFdEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
				dFdEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
				dFdEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

				dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex, locVxParamIndex + 1);
				dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex, locVxParamIndex + 2);
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 1);
			}

		}
		else //2 & 3 //px is largest
		{
			dF(locFIndex, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(locKinFitParticleType == d_DecayingParticle) //decaying particle
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 7" << endl;

				dFdXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
				dFdXi(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X();
				dFdXi(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
				dFdXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

				dFdXi(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
				dFdXi(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();
				dFdXi(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
				dFdXi(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 2);
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 1);
			}
			else
			{
				if(dDebugLevel > 30)
					cout << "DKinFitter: Calc_dF_Vertex() Section 8" << endl;

				dFdEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
				dFdEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X();
				dFdEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
				dFdEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

				dFdEta(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
				dFdEta(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();
				dFdEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
				dFdEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 2);
				dFdXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex);
				dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 1);
			}
		}
	}
}

void DKinFitter::Calc_dF_Time(size_t& locFIndex, const DKinFitParticle* locKinFitParticle, bool locUseRFTimeFlag)
{
//make sure that the particle properties of the beam used for the RF time constraint are done properly (point to correct cov matrix elements, etc.)
//for beam particles: be careful about time constraint

//charged accel: 4 constraints, 5 unknowns

//other: 3 constraints

	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
	int locCharge = locKinFitParticle->Get_Charge();
	double locShowerEnergy = locKinFitParticle->Get_ShowerEnergy();
	double locEnergy = locKinFitParticle->Get_Energy();
	double locMass = locKinFitParticle->Get_Mass();

	double locPathLength = locKinFitParticle->Get_PathLength();

	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locTParamIndex = locKinFitParticle->Get_TParamIndex();
	int locLParamIndex = locKinFitParticle->Get_LParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();
	int locCommonTParamIndex = locKinFitParticle->Get_CommonTParamIndex();

	double locTime = locUseRFTimeFlag ? dRFTime : locKinFitParticle->Get_Time();
	double locDeltaT = locKinFitParticle->Get_CommonTime() - locTime;

	double locEOverC = locEnergy/29.9792458;
	if(locUseRFTimeFlag && (locCharge != 0) && Get_IsBFieldNearBeamline()) //beam particle: charged & b-field
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Time() Section 1" << endl;

		dF(locFIndex, 0) = locDeltaT - locPathLength*locEOverC/locMomentum.Mag();

		TVector3 locMomentumTerm = locMomentum*locMass*locMass*locPathLength*(1.0/(29.9792458*locEnergy*locMomentum.Mag()*locMomentum.Mag2()));
		dFdEta(locFIndex, locPxParamIndex) = locMomentumTerm.X();
		dFdEta(locFIndex, locPxParamIndex + 1) = locMomentumTerm.Y();
		dFdEta(locFIndex, locPxParamIndex + 2) = locMomentumTerm.Z();
		dFdEta(locFIndex, locTParamIndex) = -1.0;

		dFdXi(locFIndex, locCommonTParamIndex) = 1.0;
		dFdXi(locFIndex, locLParamIndex) = -1.0*locEOverC/locMomentum.Mag();

		++locFIndex;
	}
	else if(locUseRFTimeFlag) //beam particle: neutral or no b-field
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Time() Section 2" << endl;

		double locPDotDeltaX = locDeltaX.Dot(locMomentum);
		double locPMagSq = locMomentum.Mag2();

		dF(locFIndex, 0) = locDeltaT - locEOverC*locPDotDeltaX/locPMagSq;

		TVector3 locW = locMomentum*locPDotDeltaX*(locEnergy*locEnergy + locMass*locMass)*(1.0/(29.9792458*locEnergy*locPMagSq*locPMagSq)) - locDeltaX*(locEOverC/locPMagSq);

		dFdEta(locFIndex, locPxParamIndex) = locW.X();
		dFdEta(locFIndex, locPxParamIndex + 1) = locW.Y();
		dFdEta(locFIndex, locPxParamIndex + 2) = locW.Z();

		TVector3 locVertexTerm = locMomentum*(locEOverC/locPMagSq);

		dFdEta(locFIndex, locVxParamIndex) = locVertexTerm.X();
		dFdEta(locFIndex, locVxParamIndex + 1) = locVertexTerm.Y();
		dFdEta(locFIndex, locVxParamIndex + 2) = locVertexTerm.Z();
		dFdEta(locFIndex, locTParamIndex) = -1.0;

		dFdXi(locFIndex, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex, locVxParamIndex);
		dFdXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex, locVxParamIndex + 1);
		dFdXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex, locVxParamIndex + 2);
		dFdXi(locFIndex, locCommonTParamIndex) = 1.0;

		++locFIndex;
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag()) //neutral shower
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Time() Section 3" << endl;

		dF(locFIndex, 0) = locDeltaT + locShowerEnergy*locDeltaX.Mag()/(29.9792458*locMomentum.Mag());

		dFdEta(locFIndex, locEParamIndex) = -1.0*locMass*locMass*locDeltaX.Mag()/(29.9792458*locMomentum.Mag()*locMomentum.Mag2());

		TVector3 locPositionTerm = locDeltaX.Unit()*(locEOverC/(locMomentum.Mag()));

		dFdEta(locFIndex, locVxParamIndex) = -1.0*locPositionTerm.X();
		dFdEta(locFIndex, locVxParamIndex + 1) = -1.0*locPositionTerm.Y();
		dFdEta(locFIndex, locVxParamIndex + 2) = -1.0*locPositionTerm.Z();
		dFdEta(locFIndex, locTParamIndex) = -1.0;

		dFdXi(locFIndex, locCommonVxParamIndex) = locPositionTerm.X();
		dFdXi(locFIndex, locCommonVxParamIndex + 1) = locPositionTerm.Y();
		dFdXi(locFIndex, locCommonVxParamIndex + 2) = locPositionTerm.Z();
		dFdXi(locFIndex, locCommonTParamIndex) = 1.0;

		++locFIndex;
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline()) //charged track in magnetic field
	{
		TVector3 locBField = Get_BField(locPosition);

		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locPCrossH = locMomentum.Cross(locH);
		TVector3 locPCrossHCrossH = locPCrossH.Cross(locH);
		double locPMag = locMomentum.Mag();
		double locPMagCubed = locMomentum.Mag2()*locPMag;
		double locALOverPMag = locA*locPathLength/locPMag;
		double locPDotH = locMomentum.Dot(locH);
		double locSinALOverPMag = sin(locALOverPMag);
		double locCosALOverPMag = cos(locALOverPMag);
		double locC = 1.0 - locCosALOverPMag;

		TVector3 locVertexConstraints = locPCrossHCrossH*(locSinALOverPMag/locA) + locPCrossH*(locC/locA) - (locPDotH*locPathLength/locPMag)*locH;

		dF(locFIndex, 0) = locVertexConstraints.X();
		dF(locFIndex + 1, 0) = locVertexConstraints.Y();
		dF(locFIndex + 2, 0) = locVertexConstraints.Z();
		dF(locFIndex + 3, 0) = locDeltaT - locPathLength*locEOverC/locPMag;

		//kx := k = x, etc.
		TVector3 locLPOverPCubed = locMomentum*(locPathLength/locPMagCubed);
		TVector3 locLHOverPMag = locH*(locPathLength/locPMag);
		TVector3 locYkx = locLPOverPCubed*(locPDotH*locH.X() - locPCrossHCrossH.X()*locCosALOverPMag - locPCrossH.X()*locSinALOverPMag) - locLHOverPMag*locH.X();
		TVector3 locYky = locLPOverPCubed*(locPDotH*locH.Y() - locPCrossHCrossH.Y()*locCosALOverPMag - locPCrossH.Y()*locSinALOverPMag) - locLHOverPMag*locH.Y();
		TVector3 locYkz = locLPOverPCubed*(locPDotH*locH.Z() - locPCrossHCrossH.Z()*locCosALOverPMag - locPCrossH.Z()*locSinALOverPMag) - locLHOverPMag*locH.Z();

		TVector3 locDkx = locYkx + locH*(locH.X()*locSinALOverPMag/locA);
		TVector3 locDky = locYkz + locH*(locH.Y()*locSinALOverPMag/locA);
		TVector3 locDkz = locYky + locH*(locH.Z()*locSinALOverPMag/locA);

		TVector3 locG(locYkx.X() - (1.0 - locH.X()*locH.X())*locSinALOverPMag/locA, locYky.Y() - (1.0 - locH.Y()*locH.Y())*locSinALOverPMag/locA, locYkz.Z() - (1.0 - locH.Z()*locH.Z())*locSinALOverPMag/locA);
		double locZ = locPathLength*locMass*locMass/(29.9792458*locEnergy*locPMagCubed);

		TVector3 locHlocCOverA = locH*(locC/locA);

		TVector3 locR = locDeltaX - locLHOverPMag*locPDotH + locPCrossH*(locC/locA) + locPCrossHCrossH*(locSinALOverPMag/locA);

		if(locKinFitParticleType == d_DecayingParticle) //decaying charged particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 4" << endl;

			dFdXi(locFIndex, locPxParamIndex) = locG.X();
			dFdXi(locFIndex + 1, locPxParamIndex + 1) = locG.Y();
			dFdXi(locFIndex + 2, locPxParamIndex + 2) = locG.Z();

			dFdXi(locFIndex, locPxParamIndex + 1) = locDkx.Y() + locHlocCOverA.Z();
			dFdXi(locFIndex, locPxParamIndex + 2) = locDkx.Z() - locHlocCOverA.Y();

			dFdXi(locFIndex + 1, locPxParamIndex) = locDky.X() - locHlocCOverA.Z();
			dFdXi(locFIndex + 1, locPxParamIndex + 2) = locDky.Z() + locHlocCOverA.X();

			dFdXi(locFIndex + 2, locPxParamIndex) = locDkz.X() + locHlocCOverA.Y();
			dFdXi(locFIndex + 2, locPxParamIndex + 1) = locDkz.Y() - locHlocCOverA.X();

			dFdXi(locFIndex + 3, locPxParamIndex) = locZ*locMomentum.X();
			dFdXi(locFIndex + 3, locPxParamIndex + 1) = locZ*locMomentum.Y();
			dFdXi(locFIndex + 3, locPxParamIndex + 2) = locZ*locMomentum.Z();

			dFdXi(locFIndex, locVxParamIndex) = -1.0;
			dFdXi(locFIndex + 1, locVxParamIndex + 1) = -1.0;
			dFdXi(locFIndex + 2, locVxParamIndex + 2) = -1.0;
			dFdXi(locFIndex + 3, locTParamIndex) = -1.0;

			dFdXi(locFIndex, locCommonVxParamIndex) = 1.0;
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = 1.0;
			dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) = 1.0;
			dFdXi(locFIndex + 3, locCommonTParamIndex) = 1.0;

			dFdXi(locFIndex, locLParamIndex) = locR.X();
			dFdXi(locFIndex + 1, locLParamIndex + 1) = locR.Y();
			dFdXi(locFIndex + 2, locLParamIndex + 2) = locR.Z();
			dFdXi(locFIndex + 3, locLParamIndex) = -1.0*locEOverC/locPMag;
		}
		else //detected particle or beam particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 5" << endl;

			dFdEta(locFIndex, locPxParamIndex) = locG.X();
			dFdEta(locFIndex + 1, locPxParamIndex + 1) = locG.Y();
			dFdEta(locFIndex + 2, locPxParamIndex + 2) = locG.Z();

			dFdEta(locFIndex, locPxParamIndex + 1) = locDkx.Y() + locHlocCOverA.Z();
			dFdEta(locFIndex, locPxParamIndex + 2) = locDkx.Z() - locHlocCOverA.Y();

			dFdEta(locFIndex + 1, locPxParamIndex) = locDky.X() - locHlocCOverA.Z();
			dFdEta(locFIndex + 1, locPxParamIndex + 2) = locDky.Z() + locHlocCOverA.X();

			dFdEta(locFIndex + 2, locPxParamIndex) = locDkz.X() + locHlocCOverA.Y();
			dFdEta(locFIndex + 2, locPxParamIndex + 1) = locDkz.Y() - locHlocCOverA.X();

			dFdEta(locFIndex + 3, locPxParamIndex) = locZ*locMomentum.X();
			dFdEta(locFIndex + 3, locPxParamIndex + 1) = locZ*locMomentum.Y();
			dFdEta(locFIndex + 3, locPxParamIndex + 2) = locZ*locMomentum.Z();

			dFdEta(locFIndex, locVxParamIndex) = -1.0;
			dFdEta(locFIndex + 1, locVxParamIndex + 1) = -1.0;
			dFdEta(locFIndex + 2, locVxParamIndex + 2) = -1.0;
			dFdEta(locFIndex + 3, locTParamIndex) = -1.0;

			dFdXi(locFIndex, locCommonVxParamIndex) = 1.0;
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = 1.0;
			dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) = 1.0;
			dFdXi(locFIndex + 3, locCommonTParamIndex) = 1.0;

			dFdXi(locFIndex, locLParamIndex) = locR.X();
			dFdXi(locFIndex + 1, locLParamIndex + 1) = locR.Y();
			dFdXi(locFIndex + 2, locLParamIndex + 2) = locR.Z();
			dFdXi(locFIndex + 3, locLParamIndex) = -1.0*locEOverC/locPMag;

		}
		locFIndex += 4;
	}
	else //non-accelerating charged particles, decaying neutral particles, detected neutral particles with a pre-ordained momentum
	{
		dF(locFIndex, 0) = locDeltaT - locEOverC*locDeltaX.X()/locMomentum.X();
		dF(locFIndex + 1, 0) = locDeltaT - locEOverC*locDeltaX.Y()/locMomentum.Y();
		dF(locFIndex + 2, 0) = locDeltaT - locEOverC*locDeltaX.Z()/locMomentum.Z();

		TVector3 locDeltaXOverPComponents(locDeltaX.X()/locMomentum.X(), locDeltaX.Y()/locMomentum.Y(), locDeltaX.Z()/locMomentum.Z());
		TVector3 locEOverCP(locEOverC/locMomentum.X(), locEOverC/locMomentum.Y(), locEOverC/locMomentum.Z());
		TVector3 locDeltaXOverCE = locDeltaX*(1.0/(29.9792458*locEnergy));
		TVector3 locPOverCE = locMomentum*(1.0/(29.9792458*locEnergy));

		if(locKinFitParticleType == d_DecayingParticle) //decaying particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 6" << endl;

			dFdXi(locFIndex, locPxParamIndex) = locEOverCP.X()*locDeltaXOverPComponents.X() - locDeltaXOverCE.X();
			dFdXi(locFIndex + 1, locPxParamIndex + 1) = locEOverCP.Y()*locDeltaXOverPComponents.Y() - locDeltaXOverCE.Y();
			dFdXi(locFIndex + 2, locPxParamIndex + 2) = locEOverCP.Z()*locDeltaXOverPComponents.Z() - locDeltaXOverCE.Z();

			dFdXi(locFIndex, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Y();
			dFdXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Z();

			dFdXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.X();
			dFdXi(locFIndex + 1, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.Z();

			dFdXi(locFIndex + 2, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.X();
			dFdXi(locFIndex + 2, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.Y();

			dFdXi(locFIndex, locVxParamIndex) = locEOverCP.X();
			dFdXi(locFIndex + 1, locVxParamIndex + 1) = locEOverCP.Y();
			dFdXi(locFIndex + 2, locVxParamIndex + 2) = locEOverCP.Z();

			dFdXi(locFIndex, locTParamIndex) = -1.0;
			dFdXi(locFIndex + 1, locTParamIndex) = -1.0;
			dFdXi(locFIndex + 2, locTParamIndex) = -1.0;

			dFdXi(locFIndex, locCommonVxParamIndex) = -1.0*dFdXi(locFIndex, locVxParamIndex);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdXi(locFIndex + 1, locVxParamIndex + 1);
			dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) = -1.0*dFdXi(locFIndex + 2, locVxParamIndex + 2);

			dFdXi(locFIndex, locCommonTParamIndex) = 1.0;
			dFdXi(locFIndex + 1, locCommonTParamIndex) = 1.0;
			dFdXi(locFIndex + 2, locCommonTParamIndex) = 1.0;
		}
		else //detected particle or beam particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 7" << endl;

			dFdEta(locFIndex, locPxParamIndex) = locEOverCP.X()*locDeltaXOverPComponents.X() - locDeltaXOverCE.X();
			dFdEta(locFIndex + 1, locPxParamIndex + 1) = locEOverCP.Y()*locDeltaXOverPComponents.Y() - locDeltaXOverCE.Y();
			dFdEta(locFIndex + 2, locPxParamIndex + 2) = locEOverCP.Z()*locDeltaXOverPComponents.Z() - locDeltaXOverCE.Z();

			dFdEta(locFIndex, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Y();
			dFdEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Z();

			dFdEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.X();
			dFdEta(locFIndex + 1, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.Z();

			dFdEta(locFIndex + 2, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.X();
			dFdEta(locFIndex + 2, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.Y();

			dFdEta(locFIndex, locVxParamIndex) = locEOverCP.X();
			dFdEta(locFIndex + 1, locVxParamIndex + 1) = locEOverCP.Y();
			dFdEta(locFIndex + 2, locVxParamIndex + 2) = locEOverCP.Z();

			dFdEta(locFIndex, locTParamIndex) = -1.0;
			dFdEta(locFIndex + 1, locTParamIndex) = -1.0;
			dFdEta(locFIndex + 2, locTParamIndex) = -1.0;

			dFdXi(locFIndex, locCommonVxParamIndex) = -1.0*dFdEta(locFIndex, locVxParamIndex);
			dFdXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dFdEta(locFIndex + 1, locVxParamIndex + 1);
			dFdXi(locFIndex + 2, locCommonVxParamIndex + 2) = -1.0*dFdEta(locFIndex + 2, locVxParamIndex + 2);

			dFdXi(locFIndex, locCommonTParamIndex) = 1.0;
			dFdXi(locFIndex + 1, locCommonTParamIndex) = 1.0;
			dFdXi(locFIndex + 2, locCommonTParamIndex) = 1.0;

		}
		locFIndex += 3;
	}
}

void DKinFitter::Calc_Pulls(void)
{
	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	int locParamIndex;
	double locDenominator;

	locParamIndex = 0;
	dPulls.clear();
	map<DKinFitPullType, double> locParticlePulls;
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticle = dKinFitParticles[loc_i];
		locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;
		locParticlePulls.clear();

		if(dDebugLevel >= 50)
		{
			cout << "pulls: q, mass = " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
			cout << "e, px, vx, t param indices = " << dKinFitParticles[loc_i]->Get_EParamIndex() << ", " << dKinFitParticles[loc_i]->Get_PxParamIndex() << ", " << dKinFitParticles[loc_i]->Get_VxParamIndex() << ", " << dKinFitParticles[loc_i]->Get_TParamIndex() << endl;
		}

		locParamIndex = locKinFitParticle->Get_EParamIndex();
		if(locParamIndex >= 0) //E
		{
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_EPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		locParamIndex = locKinFitParticle->Get_PxParamIndex();
		if(locParamIndex >= 0) //px, py, pz
		{
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_PxPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
			++locParamIndex;
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_PyPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
			++locParamIndex;
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_PzPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		locParamIndex = locKinFitParticle->Get_VxParamIndex();
		if(locParamIndex >= 0) //vx, vy, vz
		{
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_XxPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
			++locParamIndex;
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_XyPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
			++locParamIndex;
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_XzPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		locParamIndex = locKinFitParticle->Get_TParamIndex();
		if(locParamIndex >= 0) //T
		{
			locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_TPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		if(!locParticlePulls.empty())
			dPulls[locKinFitParticle] = locParticlePulls;
	}

	//RF BUNCH
	locParamIndex = dRFTimeParamIndex;
	if(locParamIndex >= 0)
	{
		locParticlePulls.clear();
		locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
		locParticlePulls[d_TPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		dPulls[NULL] = locParticlePulls;
	}

	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dEpsilon: " << endl;
		Print_Matrix(dEpsilon);
		cout << "DKinFitter: dVY: " << endl;
		Print_Matrix(dVY);
		cout << "DKinFitter: dVEta: " << endl;
		Print_Matrix(*dVEta);
		cout << "DKinFitter: Pulls: " << endl;
		map<const DKinFitParticle*, map<DKinFitPullType, double> >::iterator locIterator;
		map<DKinFitPullType, double>::iterator locIterator2;
		for(locIterator = dPulls.begin(); locIterator != dPulls.end(); ++locIterator)
		{
			map<DKinFitPullType, double>& locTempParticlePulls = locIterator->second;
			const DKinFitParticle* locTempKinFitParticle = locIterator->first;
			TVector3 locMomentum = locTempKinFitParticle->Get_Momentum();
			cout << "particle q, p3 = " << locTempKinFitParticle->Get_Charge() << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ":" << endl;
			for(size_t loc_i = 0; loc_i < 8; ++loc_i)
			{
				if(locTempParticlePulls.find((DKinFitPullType)loc_i) != locTempParticlePulls.end())
					cout << locTempParticlePulls[(DKinFitPullType)loc_i] << ", ";
			}
			cout << endl;
		}
	}
}

void DKinFitter::Set_FinalTrackInfo(void)
{
	TVector3 locMomentum;
	TLorentzVector locSpacetimeVertex;
	DKinFitParticleType locKinFitParticleType;
	int locPxParamIndex, locVxParamIndex, locTParamIndex, locEParamIndex;
	int locCovMatrixEParamIndex, locCovMatrixPxParamIndex, locCovMatrixVxParamIndex, locCovMatrixTParamIndex;
	double locUncertaintyRatio;

	// first update the covariance matrices of each particle with the fit results (prior to any propagation)
		//correlations between fit and unfit measured parameters: 
			//assume that the correlations remain unchanged: the changes in the uncertainties should be small wrst the uncertainties themselves
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = dKinFitParticles[loc_i]->Get_KinFitParticleType();
		if(locKinFitParticleType == d_TargetParticle)
			continue;

		bool locReconstructedParticleFlag = ((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle));
		TMatrixDSym& locKinFitMatrix = locReconstructedParticleFlag ? *dVXi : *dVEta;
		TMatrixDSym& locCovarianceMatrix = locReconstructedParticleFlag ? *(Get_MatrixDSymResource()) : *(dKinFitParticles[loc_i]->dCovarianceMatrix);
		if(locReconstructedParticleFlag)
			locCovarianceMatrix.ResizeTo(7, 7);

		locPxParamIndex = dKinFitParticles[loc_i]->Get_PxParamIndex();
		locVxParamIndex = dKinFitParticles[loc_i]->Get_VxParamIndex();
		locTParamIndex = dKinFitParticles[loc_i]->Get_TParamIndex();
		locEParamIndex = dKinFitParticles[loc_i]->Get_EParamIndex();

		locCovMatrixEParamIndex = dKinFitParticles[loc_i]->Get_CovMatrixEParamIndex();
		locCovMatrixPxParamIndex = dKinFitParticles[loc_i]->Get_CovMatrixPxParamIndex();
		locCovMatrixVxParamIndex = dKinFitParticles[loc_i]->Get_CovMatrixVxParamIndex();
		locCovMatrixTParamIndex = dKinFitParticles[loc_i]->Get_CovMatrixTParamIndex();

		if(dDebugLevel >= 50)
		{
			cout << "SETTING FINAL TRACK INFO: q, mass = " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
			cout << "E, px, vx, t param indices = " << locEParamIndex << ", " << locPxParamIndex << ", " << locVxParamIndex << ", " << locTParamIndex << endl;
			cout << "E, px, vx, t cov matrix indices = " << locCovMatrixEParamIndex << ", " << locCovMatrixPxParamIndex << ", " << locCovMatrixVxParamIndex << ", " << locCovMatrixTParamIndex << endl;
			cout << "sizes = " << locCovarianceMatrix.GetNrows() << ", " << locCovarianceMatrix.GetNcols() << ", " << locKinFitMatrix.GetNrows() << ", " << locKinFitMatrix.GetNcols() << endl;
		}

		//set localized terms first
		if(locEParamIndex >= 0)
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixEParamIndex) = locKinFitMatrix(locEParamIndex, locEParamIndex);
		if(locPxParamIndex >= 0)
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
					locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixPxParamIndex + loc_k) = locKinFitMatrix(loc_j + locPxParamIndex, loc_k + locPxParamIndex);
			}
		}
		if(locVxParamIndex >= 0)
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
					locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixVxParamIndex + loc_k) = locKinFitMatrix(loc_j + locVxParamIndex, loc_k + locVxParamIndex);
			}
		}
		if(locTParamIndex >= 0)
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixTParamIndex) = locKinFitMatrix(locTParamIndex, locTParamIndex);

		//cross terms
		if((locEParamIndex >= 0) && (locVxParamIndex >= 0)) //both included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixEParamIndex + 0, locCovMatrixVxParamIndex + loc_j) = locKinFitMatrix(locEParamIndex + 0, locVxParamIndex + loc_j);
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixEParamIndex + 0) = locKinFitMatrix(locVxParamIndex + loc_j, locEParamIndex + 0);
			}
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex >= 0) && (locVxParamIndex < 0)) //only E included in the fit
		{
			locUncertaintyRatio = sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/sqrt(dVY(locEParamIndex, locEParamIndex));
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixEParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixEParamIndex + 0) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex < 0) && (locVxParamIndex >= 0) && (locCovMatrixEParamIndex >= 0)) //only V included in the fit //E may not be in the covariance matrix!!
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locUncertaintyRatio = sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				locCovarianceMatrix(locCovMatrixEParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixEParamIndex + 0) *= locUncertaintyRatio;
			}
		}

		if((locEParamIndex >= 0) && (locTParamIndex >= 0)) //both included in the fit
		{
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) = locKinFitMatrix(locEParamIndex, locTParamIndex);
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) = locKinFitMatrix(locTParamIndex, locEParamIndex);
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex >= 0) && (locTParamIndex < 0)) //only E included in the fit
		{
			locUncertaintyRatio = sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/sqrt(dVY(locEParamIndex, locEParamIndex));
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) *= locUncertaintyRatio;
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) *= locUncertaintyRatio;
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex < 0) && (locTParamIndex >= 0) && (locCovMatrixEParamIndex >= 0)) //only T included in the fit //E may not be in the covariance matrix!!
		{
			locUncertaintyRatio = sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/sqrt(dVY(locTParamIndex, locTParamIndex));
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) *= locUncertaintyRatio;
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) *= locUncertaintyRatio;
		}

		if((locPxParamIndex >= 0) && (locVxParamIndex >= 0)) //both included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				{
					locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixVxParamIndex + loc_k) = locKinFitMatrix(locPxParamIndex + loc_j, locVxParamIndex + loc_k);
					locCovarianceMatrix(locCovMatrixVxParamIndex + loc_k, locCovMatrixPxParamIndex + loc_j) = locKinFitMatrix(locVxParamIndex + loc_k, locPxParamIndex + loc_j);
				}
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex >= 0) && (locVxParamIndex < 0)) //only P included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locUncertaintyRatio = sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/sqrt(dVY(locPxParamIndex + loc_j, locPxParamIndex + loc_j));
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				{
					locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixVxParamIndex + loc_k) *= locUncertaintyRatio;
					locCovarianceMatrix(locCovMatrixVxParamIndex + loc_k, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
				}
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex < 0) && (locVxParamIndex >= 0)) //only V included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locUncertaintyRatio = sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				{
					locCovarianceMatrix(locCovMatrixPxParamIndex + loc_k, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
					locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixPxParamIndex + loc_k) *= locUncertaintyRatio;
				}
			}
		}

		if((locPxParamIndex >= 0) && (locTParamIndex >= 0)) //both included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) = locKinFitMatrix(locPxParamIndex + loc_j, locTParamIndex + 0);
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) = locKinFitMatrix(locTParamIndex + 0, locPxParamIndex + loc_j);
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex >= 0) && (locTParamIndex < 0)) //only P included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locUncertaintyRatio = sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/sqrt(dVY(locPxParamIndex + loc_j, locPxParamIndex + loc_j));
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			locUncertaintyRatio = sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/sqrt(dVY(locTParamIndex, locTParamIndex));
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}

		if((locVxParamIndex >= 0) && (locTParamIndex >= 0)) //both included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) = locKinFitMatrix(locVxParamIndex + loc_j, locTParamIndex + 0);
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) = locKinFitMatrix(locTParamIndex + 0, locVxParamIndex + loc_j);
			}
		}
		else if(!locReconstructedParticleFlag && (locVxParamIndex >= 0) && (locTParamIndex < 0)) //only V included in the fit
		{
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locUncertaintyRatio = sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locVxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			locUncertaintyRatio = sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/sqrt(dVY(locTParamIndex, locTParamIndex));
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
	}

	// propagate the track parameters
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = dKinFitParticles[loc_i]->Get_KinFitParticleType();

		if(!dKinFitParticles[loc_i]->Get_IsInVertexOrSpacetimeFitFlag())
			continue; // no distance over which to propagate

		if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			continue; // particle properties already defined at the fit vertex

		pair<double, double> locPathLengthPair;
		TMatrixDSym& locCovarianceMatrix = *(dKinFitParticles[loc_i]->dCovarianceMatrix);

		if(!Propagate_TrackInfoToCommonVertex(dKinFitParticles[loc_i], dVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, locCovarianceMatrix))
			continue; // info not propagated

		if(dDebugLevel >= 50)
		{
			cout << "PROPAGATED FINAL TRACK INFO: q, mass = " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
			cout << "p_xyz, v_xyzt = " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;
			cout << "common v_xyzt = " << dKinFitParticles[loc_i]->Get_CommonVertex().X() << ", " << dKinFitParticles[loc_i]->Get_CommonVertex().Y() << ", " << dKinFitParticles[loc_i]->Get_CommonVertex().Z() << ", " << dKinFitParticles[loc_i]->Get_CommonTime() << endl;
			cout << "path length & uncert = " << locPathLengthPair.first << ", " << locPathLengthPair.second << endl;
			cout << "sizes = " << locCovarianceMatrix.GetNrows() << ", " << locCovarianceMatrix.GetNcols() << endl;
		}

		//no need to set the covariance matrix: already updated (passed in a reference to it)
		dKinFitParticles[loc_i]->Set_Momentum(locMomentum);
		dKinFitParticles[loc_i]->Set_Position(locSpacetimeVertex.Vect());
		dKinFitParticles[loc_i]->Set_Time(locSpacetimeVertex.T());
		dKinFitParticles[loc_i]->Set_PathLength(locPathLengthPair.first);
		dKinFitParticles[loc_i]->Set_PathLengthUncertainty(locPathLengthPair.second);
	}

	//calculate the path length of decaying particles involved in 2 vertex fits
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = dKinFitParticles[loc_i]->Get_KinFitParticleType();

		if((locKinFitParticleType != d_DecayingParticle) || (dKinFitParticles[loc_i]->Get_NumVertexFits() != 2))
			continue;

		pair<double, double> locPathLengthPair;
		const TMatrixDSym& locCovarianceMatrix = *(dKinFitParticles[loc_i]->Get_CovarianceMatrix());
		if(Calc_PathLength(dKinFitParticles[loc_i], dVXi, locCovarianceMatrix, locPathLengthPair))
		{
			dKinFitParticles[loc_i]->Set_PathLength(locPathLengthPair.first);
			dKinFitParticles[loc_i]->Set_PathLengthUncertainty(locPathLengthPair.second);
		}
	}
}

// propagates the track info to the fit common vertex
	//returns false if nothing changed (info not propagated: e.g. missing particle), else returns true
//assumes: that between the two points on the track (input measured point & kinfit common point): the b-field is constant and there is no eloss or multiple scattering 
// this function acts in the following way on each particle, IN THIS ORDER OF EXECUTION:
	// if a common vertex was not fit, then the current results are returned
	// if this is either a target particle, a missing particle, or a decaying particle that is present in only one vertex fit: 
		// the current info is defined to be at the common vertex: just return it
	// if this is a decaying particle that is present in two vertex fits: the returned info is the info at the "other" vertex
		// "other" is the opposite of where the particle info is defined to be: 
		// e.g. if dDecayingParticleAtProductionVertexFlag == true, then the "other" vertex is the decay vertex
	// for the remaining types (including decaying particles):
		// the vertex is set to the common vertex
		// if the time was fit, then the time is set to the common time, else it is propagated to the common vertex
		// if the path length was fit then the fit result is used; otherwise it is propagated
		// momentum:
			// if this is a neutral shower: momentum is redefined by the new vertex //already done by Update_ParticleParams()
			// if this is a neutral particle (either detected, beam, or decaying) or a charged particle without a b-field: the momentum is unchanged
			// if this is a charged particle in a b-field (either beam, detected, or decaying): the momentum is propagated to the common vertex
	//propagating the covariance matrix:
		//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
		//add common time to matrix: 11x11 or 9x9 (neutral shower): if kinfit just add in (no correlations to meas, corr to common v3), else transform
		//transform to 7x7: common v3 & common t are just copied to the measured spots
//the output covariance matrix is 7x7, even if the particle represents a neutral shower (5x5)
//the only thing this function uses that depends on the state of the fitter itself is the magnetic field map
bool DKinFitter::Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixDSym& locCovarianceMatrix) const
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if(!locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag())
		return false; // no distance over which to propagate

	if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
		return false; // particle properties already defined at the fit vertex

	if((locKinFitParticleType == d_DecayingParticle) && (locKinFitParticle->Get_NumVertexFits() != 2))
		return false; // particle properties already defined at the fit vertex

	locCovarianceMatrix = *(locKinFitParticle->Get_CovarianceMatrix());
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locBField = Get_BField(locPosition);
	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
	double locCommonTime;

	// covariance matrix
	int locCovMatrixEParamIndex = locKinFitParticle->Get_CovMatrixEParamIndex();
	int locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
	int locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
	int locCovMatrixTParamIndex = locKinFitParticle->Get_CovMatrixTParamIndex();
	int locCommonVxParamIndex_TempMatrix, locCommonTParamIndex_TempMatrix;

	DKinFitConstraint_VertexBase* locKinFitConstraint_VertexBase = locKinFitParticle->Get_CommonVertexAndOrTimeConstraint();
	int locCommonVxParamIndex = locKinFitConstraint_VertexBase->Get_VxParamIndex();
	int locCommonTParamIndex = -1;
	if(locKinFitParticle->Get_IsInSpacetimeFitFlag()) //vertex & time were fit
		locCommonTParamIndex = (static_cast<DKinFitConstraint_Spacetime*>(locKinFitConstraint_VertexBase))->Get_TParamIndex();

	//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
	locCommonVxParamIndex_TempMatrix = locCovarianceMatrix.GetNcols();
	locCovarianceMatrix.ResizeTo(locCommonVxParamIndex_TempMatrix + 3, locCommonVxParamIndex_TempMatrix + 3);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locCovarianceMatrix(loc_i + locCommonVxParamIndex_TempMatrix, loc_j + locCommonVxParamIndex_TempMatrix) = (*locVXi)(locCommonVxParamIndex + loc_i, locCommonVxParamIndex + loc_j);
	}
	// done: no correlations between common vertex and measured params!!!
	locSpacetimeVertex.SetVect(locCommonVertex);

	//add common time to matrix: 11x11 or 9x9 (neutral shower): if kinfit just add in (no correlations to meas, corr to common v3), else transform
		//note that if the common time is not kinfit, then the true uncertainty is overestimated: cannot be obtained without a kinematic fit (3 equations (xyz), one unknown (time))
	locCommonTParamIndex_TempMatrix = locCovarianceMatrix.GetNcols();
	if(locKinFitParticle->Get_IsInSpacetimeFitFlag()) //spacetime was fit
	{
		locCommonTime = locKinFitParticle->Get_CommonTime();
		locCovarianceMatrix.ResizeTo(locCovarianceMatrix.GetNcols() + 1, locCovarianceMatrix.GetNcols() + 1);
		locCovarianceMatrix(locCommonTParamIndex_TempMatrix, locCommonTParamIndex_TempMatrix) = (*locVXi)(locCommonTParamIndex, locCommonTParamIndex);
		for(size_t loc_i = 0; loc_i < 3; ++loc_i) //correlations to common v3
		{
			locCovarianceMatrix(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + loc_i) = (*locVXi)(locCommonTParamIndex, locCommonVxParamIndex + loc_i);
			locCovarianceMatrix(locCommonVxParamIndex_TempMatrix + loc_i, locCommonTParamIndex_TempMatrix) = (*locVXi)(locCommonVxParamIndex + loc_i, locCommonTParamIndex);
		}
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locP4.Vect().Dot(locH);
		locCommonTime = locKinFitParticle->Get_Time() + locDeltaXDotH*locP4.E()/(29.9792458*locPDotH);

		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix.GetNcols() + 1, locCovarianceMatrix.GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix.GetNcols(); ++loc_i)
			locTransformationMatrix_CommonTime(loc_i, loc_i) = 1.0; //other params are unchanged

		TVector3 locDCommonTimeDP3 = (locDeltaXDotH/(29.9792458*locPDotH)) * ((1.0/locP4.E())*locP4.Vect() - (locP4.E()/locPDotH)*locH);
		TVector3 locDCommonTimeDCommonVertex = (locP4.E()/(29.9792458*locPDotH))*locH;
		TVector3 locDCommonTimeDPosition = -1.0*locDCommonTimeDCommonVertex;

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex) = locDCommonTimeDP3.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex + 1) = locDCommonTimeDP3.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex + 2) = locDCommonTimeDP3.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex) = locDCommonTimeDPosition.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 1) = locDCommonTimeDPosition.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 2) = locDCommonTimeDPosition.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix) = locDCommonTimeDCommonVertex.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 1) = locDCommonTimeDCommonVertex.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 2) = locDCommonTimeDCommonVertex.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixTParamIndex) = 1.0;

		locCovarianceMatrix.Similarity(locTransformationMatrix_CommonTime);
	}
	else if(!locNeutralShowerFlag) //non-accelerating, non-shower
	{
		double locDeltaXDotP = locDeltaX.Dot(locP4.Vect());
		locCommonTime = locKinFitParticle->Get_Time() + locDeltaXDotP*locP4.E()/(29.9792458*locP4.Vect().Mag2());

		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix.GetNcols() + 1, locCovarianceMatrix.GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix.GetNcols(); ++loc_i)
			locTransformationMatrix_CommonTime(loc_i, loc_i) = 1.0; //other params are unchanged

		TVector3 locDCommonTimeDP3 = (1.0/(29.9792458*locP4.Vect().Mag2())) * (locP4.E()*locDeltaX + locDeltaXDotP*(1.0/locP4.E() - 2.0*locP4.E()/locP4.Vect().Mag2())*locP4.Vect());
		TVector3 locDCommonTimeDCommonVertex = (locP4.E()/(29.9792458*locP4.Vect().Mag2()))*locP4.Vect();
		TVector3 locDCommonTimeDPosition = -1.0*locDCommonTimeDCommonVertex;

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex) = locDCommonTimeDP3.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex + 1) = locDCommonTimeDP3.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixPxParamIndex + 2) = locDCommonTimeDP3.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex) = locDCommonTimeDPosition.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 1) = locDCommonTimeDPosition.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 2) = locDCommonTimeDPosition.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix) = locDCommonTimeDCommonVertex.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 1) = locDCommonTimeDCommonVertex.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 2) = locDCommonTimeDCommonVertex.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixTParamIndex) = 1.0;

		locCovarianceMatrix.Similarity(locTransformationMatrix_CommonTime);
	}
	else //neutral shower
	{
		locCommonTime = locKinFitParticle->Get_Time() - locDeltaX.Mag()*locP4.E()/(29.9792458*locP4.P());

		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix.GetNcols() + 1, locCovarianceMatrix.GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix.GetNcols(); ++loc_i)
			locTransformationMatrix_CommonTime(loc_i, loc_i) = 1.0; //other params are unchanged

		double locDCommonTimeDEnergy = locDeltaX.Mag()*locP4.M2()/(29.9792458*locP4.P()*locP4.Vect().Mag2());
		TVector3 locDCommonTimeDPosition = (locP4.E()/(29.9792458*locP4.P()*locDeltaX.Mag()))*locDeltaX;
		TVector3 locDCommonTimeDCommonVertex = -1.0*locDCommonTimeDPosition;

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixEParamIndex) = locDCommonTimeDEnergy;

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex) = locDCommonTimeDPosition.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 1) = locDCommonTimeDPosition.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixVxParamIndex + 2) = locDCommonTimeDPosition.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix) = locDCommonTimeDCommonVertex.X();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 1) = locDCommonTimeDCommonVertex.Y();
		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + 2) = locDCommonTimeDCommonVertex.Z();

		locTransformationMatrix_CommonTime(locCommonTParamIndex_TempMatrix, locCovMatrixTParamIndex) = 1.0;

		locCovarianceMatrix.Similarity(locTransformationMatrix_CommonTime);
	}
	locSpacetimeVertex.SetT(locCommonTime);

	//transform to 7x7: common v3 & common t are just copied to the measured spots; p is propagated if in bfield, else is copied
	TMatrixD locTransformationMatrix_Propagation(7, locCovarianceMatrix.GetNcols());
	//p3
	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //charged & in b-field
	{
		locMomentum = locP4.Vect() - locDeltaX.Cross(locA*locH);

		locTransformationMatrix_Propagation(0, locCovMatrixPxParamIndex) = 1.0;
		locTransformationMatrix_Propagation(1, locCovMatrixPxParamIndex + 1) = 1.0;
		locTransformationMatrix_Propagation(2, locCovMatrixPxParamIndex + 2) = 1.0;

		locTransformationMatrix_Propagation(0, locCovMatrixVxParamIndex + 1) = -1.0*locA*locH.Z();
		locTransformationMatrix_Propagation(0, locCovMatrixVxParamIndex + 2) = locA*locH.Y();

		locTransformationMatrix_Propagation(1, locCovMatrixVxParamIndex) = locA*locH.Z();
		locTransformationMatrix_Propagation(1, locCovMatrixVxParamIndex + 2) = -1.0*locA*locH.X();

		locTransformationMatrix_Propagation(2, locCovMatrixVxParamIndex) = -1.0*locA*locH.Y();
		locTransformationMatrix_Propagation(2, locCovMatrixVxParamIndex + 1) = locA*locH.X();

		locTransformationMatrix_Propagation(0, locCommonVxParamIndex_TempMatrix + 1) = locA*locH.Z();
		locTransformationMatrix_Propagation(0, locCommonVxParamIndex_TempMatrix + 2) = -1.0*locA*locH.Y();

		locTransformationMatrix_Propagation(1, locCommonVxParamIndex_TempMatrix) = -1.0*locA*locH.Z();
		locTransformationMatrix_Propagation(1, locCommonVxParamIndex_TempMatrix + 2) = locA*locH.X();

		locTransformationMatrix_Propagation(2, locCommonVxParamIndex_TempMatrix) = locA*locH.Y();
		locTransformationMatrix_Propagation(2, locCommonVxParamIndex_TempMatrix + 1) = -1.0*locA*locH.X();
	}
	else //constant: either neutral or no b-field
	{
		locMomentum = locP4.Vect();
		for(unsigned int loc_i = 0; loc_i < 3; ++loc_i)
			locTransformationMatrix_Propagation(loc_i, locCovMatrixPxParamIndex + loc_i) = 1.0;
	}
	//v3
	for(unsigned int loc_i = 0; loc_i < 3; ++loc_i)
		locTransformationMatrix_Propagation(3 + loc_i, locCommonVxParamIndex_TempMatrix + loc_i) = 1.0;
	//t
	locTransformationMatrix_Propagation(6, locCommonTParamIndex_TempMatrix) = 1.0;
	//transform!!
	locCovarianceMatrix.Similarity(locTransformationMatrix_Propagation); //FINALLY!!!

	return Calc_PathLength(locKinFitParticle, locVXi, locCovarianceMatrix, locPathLengthPair);
}

bool DKinFitter::Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixDSym& locCovarianceMatrix, pair<double, double>& locPathLengthPair) const
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if(!locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag())
		return false; // no distance over which to propagate

	if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
		return false; // particle properties already defined at the fit vertex

	if((locKinFitParticleType == d_DecayingParticle) && (locKinFitParticle->Get_NumVertexFits() != 2))
		return false; // particle properties already defined at the fit vertex

	double locPathLength, locPathLengthUncertainty;
	int locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
	int locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
	DKinFitConstraint_VertexBase* locKinFitConstraint_VertexBase = locKinFitParticle->Get_CommonVertexAndOrTimeConstraint();
	int locCommonVxParamIndex = locKinFitConstraint_VertexBase->Get_VxParamIndex();

	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locBField = Get_BField(locPosition);
	TVector3 locH = locBField.Unit();

	int locLParamIndex = locKinFitParticle->Get_LParamIndex();

	//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
	TMatrixDSym locTempMatrix(locCovarianceMatrix);
	int locCommonVxParamIndex_TempMatrix = locTempMatrix.GetNcols();
	locTempMatrix.ResizeTo(locCommonVxParamIndex_TempMatrix + 3, locCommonVxParamIndex_TempMatrix + 3);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			locTempMatrix(loc_i + locCommonVxParamIndex_TempMatrix, loc_j + locCommonVxParamIndex_TempMatrix) = (*locVXi)(locCommonVxParamIndex + loc_i, locCommonVxParamIndex + loc_j);
	}

	//find path length & its uncertainty
	if(locKinFitParticle->Get_IsInSpacetimeFitFlag() && (locCharge != 0) && Get_IsBFieldNearBeamline()) //path length was fit
	{
		locPathLength = locKinFitParticle->Get_PathLength();
		locPathLengthUncertainty = (*locVXi)(locLParamIndex, locLParamIndex);
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locP4.Vect().Dot(locH);
		double locPMag = locP4.P();
		locPathLength = locDeltaXDotH*locPMag/locPDotH; //cos(theta_helix) = p.dot(h)/|p| = x.dot(h)/l (l = path length)

		TMatrixD locTransformationMatrix_PathLength(1, locTempMatrix.GetNcols());

		TVector3 locDPathDP3 = (locDeltaXDotH/locPDotH)*((1.0/locPMag)*locP4.Vect() - (locPMag/locPDotH)*locH);
		TVector3 locDPathDCommon = (locPMag/locPDotH)*locH;
		TVector3 locDPathDPosition = -1.0*locDPathDCommon;

		locTransformationMatrix_PathLength(0, locCovMatrixPxParamIndex) = locDPathDP3.X();
		locTransformationMatrix_PathLength(0, locCovMatrixPxParamIndex + 1) = locDPathDP3.Y();
		locTransformationMatrix_PathLength(0, locCovMatrixPxParamIndex + 2) = locDPathDP3.Z();

		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex) = locDPathDPosition.X();
		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex + 1) = locDPathDPosition.Y();
		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex + 2) = locDPathDPosition.Z();

		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix) = locDPathDCommon.X();
		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix + 1) = locDPathDCommon.Y();
		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix + 2) = locDPathDCommon.Z();

		locTempMatrix.Similarity(locTransformationMatrix_PathLength);
		locPathLengthUncertainty = sqrt(locTempMatrix(0, 0));
	}
	else // non-accelerating
	{
		locPathLength = locDeltaX.Mag();

		TMatrixD locTransformationMatrix_PathLength(1, locTempMatrix.GetNcols());

		TVector3 locDPathDCommon = locDeltaX.Unit();
		TVector3 locDPathDPosition = -1.0*locDPathDCommon;

		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex) = locDPathDPosition.X();
		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex + 1) = locDPathDPosition.Y();
		locTransformationMatrix_PathLength(0, locCovMatrixVxParamIndex + 2) = locDPathDPosition.Z();

		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix) = locDPathDCommon.X();
		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix + 1) = locDPathDCommon.Y();
		locTransformationMatrix_PathLength(0, locCommonVxParamIndex_TempMatrix + 2) = locDPathDCommon.Z();

		locTempMatrix.Similarity(locTransformationMatrix_PathLength);
		locPathLengthUncertainty = sqrt(locTempMatrix(0, 0));
	}
	locPathLengthPair.first = locPathLength;
	locPathLengthPair.second = locPathLengthUncertainty;

	return true;
}

