#include "DKinFitter.h"

//THINGS TO DO:
	//double check: spacetime eqs & partial derivatives
	//double check: Spacetime being set for noconstrain particles
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
	dLinkVerticesFlag = true;
	dDebugLevel = 0;

	dMaxKinFitParticlePoolSize = 100;
	dMaxKinFitConstraintVertexPoolSize = 25;
	dMaxKinFitConstraintSpacetimePoolSize = 25;
	dMaxKinFitConstraintP4PoolSize = 25;
	dMaxMatrixDSymPoolSize = 100;

	dMaxNumIterations = 20;

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

	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedInitialParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedInitialParticles[loc_i]->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one vertex constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two vertex constraints.  Constraint not added." << endl;
			return NULL;
		}
	}
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedFinalParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedFinalParticles[loc_i]->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one vertex constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two vertex constraints.  Constraint not added." << endl;
			return NULL;
		}
	}

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

	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedInitialParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedInitialParticles[loc_i]->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one spacetime constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two spacetime constraints.  Constraint not added." << endl;
			return NULL;
		}
	}
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedFinalParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedFinalParticles[loc_i]->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one spacetime constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two spacetime constraints.  Constraint not added." << endl;
			return NULL;
		}
	}

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

	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedInitialParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedInitialParticles[loc_i]->Get_NumP4Constraints();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one P4 constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two P4 constraints.  Constraint not added." << endl;
			return NULL;
		}
	}
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
	{
		DKinFitParticleType locKinFitParticleType = locClonedFinalParticles[loc_i]->Get_KinFitParticleType();
		size_t locNumP4Constraints = locClonedFinalParticles[loc_i]->Get_NumP4Constraints();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one P4 constraint.  Constraint not added." << endl;
			return NULL;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two P4 constraints.  Constraint not added." << endl;
			return NULL;
		}
	}

	//create the constraint and set its members
	DKinFitConstraint_P4* locKinFitConstraint = Get_KinFitConstraintP4Resource();
	locKinFitConstraint->Set_InitialParticles(locClonedInitialParticles);
	locKinFitConstraint->Set_FinalParticles(locClonedFinalParticles);

	//mark constraint in particles
	for(size_t loc_i = 0; loc_i < locClonedInitialParticles.size(); ++loc_i)
		locClonedInitialParticles[loc_i]->Add_P4Constraint(locKinFitConstraint);
	for(size_t loc_i = 0; loc_i < locClonedFinalParticles.size(); ++loc_i)
		locClonedFinalParticles[loc_i]->Add_P4Constraint(locKinFitConstraint);

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

bool DKinFitter::Fit_Reaction(void)
{
	if(!Resolve_Constraints())
		return false;

	Set_MatrixSizes();
	Resize_Matrices();
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
	TMatrixD locR(dNumF, 1);
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
			cout << "DKinFitter: dF_dXi: " << endl;
			Print_Matrix(dF_dXi);
			cout << "DKinFitter: dF_dEta: " << endl;
			Print_Matrix(dF_dEta);
		}

		locR = dF + dF_dEta*(dY - dEta);

		if(!Calc_dS())
		{
			if(dDebugLevel > 10)
				cout << "DKinFitter: Failed S-matrix inversion. Returning false." << endl;
			return false; // matrix is not invertible
		}

		if(dNumXi > 0)
		{
			if(!Calc_dU())
			{
				if(dDebugLevel > 10)
					cout << "DKinFitter: Failed VXi-matrix inversion. Returning false." << endl;
				return false; // matrix is not invertible
			}

			TMatrixD locDeltaXi(dNumXi, 1);
			locDeltaXi = -1.0*dU*dF_dXi_T*dS_Inverse*locR;

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

			dLambda = dS_Inverse*(locR + dF_dXi*locDeltaXi);
		}
		else
		{
			dLambda = dS_Inverse*locR;
		}

		dLambda_T.Transpose(dLambda);

		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dLambda: " << endl;
			Print_Matrix(dLambda);
		}

		dEta = dY - dVY*dF_dEta_T*dLambda;
		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dEta: " << endl;
			Print_Matrix(dEta);
		}

		TMatrixDSym locTempMatrix5 = dS;
		dChiSq = (locTempMatrix5.SimilarityT(dLambda) + 2.0*dLambda_T*dF)(0, 0);

		if(dDebugLevel > 20)
			cout << "DKinFitter: dChiSq = " << dChiSq << endl;

		Update_ParticleParams(); //input eta & xi info into particle objects
		if(dDebugLevel > 20)
		{
			for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
				Print_ParticleParams(dKinFitParticles[loc_i]);
		}
	}
	while((fabs(dChiSq - locPreviousChiSq) > 0.001) || (dChiSq < 0.0));

	// dVXi
	if(dNumXi > 0)
		*dVXi = dU;

	// dVEta
	TMatrixDSym locG = dS_Inverse;
	locG.SimilarityT(dF_dEta);
	if(dNumXi > 0)
	{
		TMatrixD locH = dF_dEta_T*dS_Inverse*dF_dXi;
		TMatrixDSym locTempMatrix11 = *dVXi;
		*dVEta = dVY - (locG - locTempMatrix11.Similarity(locH)).Similarity(dVY);
	}
	else
	{
		*dVEta = dVY - locG.Similarity(dVY); //destroys locG, but it's not needed anymore
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
	locTempMatrix.Similarity(dF_dEta);
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

bool DKinFitter::Calc_dU(void)
{
	TMatrixDSym locTempMatrix = dS_Inverse;
	locTempMatrix.SimilarityT(dF_dXi);
	dU_Inverse = locTempMatrix;
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dU_Inverse: " << endl;
		Print_Matrix(dU_Inverse);
		cout << "determinant magnitude = " << fabs(dU_Inverse.Determinant()) << endl;
	}
	TDecompLU locDecompLU_VXiInv(dU_Inverse);
	//check to make sure that the matrix is decomposable and has a non-zero determinant
	if((!locDecompLU_VXiInv.Decompose()) || (fabs(dU_Inverse.Determinant()) < 1.0E-300))
	{
		if(dDebugLevel > 10)
			cout << "DKinFitter: dU_Inverse not invertible.  Returning false." << endl;
		return false; // matrix is not invertible
	}
	dU = dU_Inverse;
	dU.Invert();
	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dU: " << endl;
		Print_Matrix(dU);
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
	cout << "DKinFitter: Particle Type Enum: " << locKinFitParticle->Get_KinFitParticleType() << endl;
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

	if(!Resolve_P4Constraints())
		return false;

	if(!Resolve_DecayingParticleSpacetimeLinks())
		return false;

	return true;
}

bool DKinFitter::Resolve_P4Constraints(void)
{
	//snag unconstrained particles, and check whether decaying particle used in vertex constraints correctly
	deque<DKinFitParticle*> locUnconstrainedParticles;
	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
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
			if(dDebugLevel > 5)
			{
				if((locConstrainableParticles.back()).first->Get_DecayingParticleAtProductionVertexFlag())
					cout << "particle is defined at it's production vertex (possibly constrained to its decay vertex)" << endl;
				else
					cout << "particle is defined at it's decay vertex (possibly constrained to its production vertex)" << endl;
			}
			locConstrainedParticles.push_back(locConstrainableParticles.back().first);
			locConstrainableParticles.pop_back();
		}
	}

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
		if(locKinFitConstraint_P4->dConstrainedP4Particle != NULL)
			continue; //constraint already set
		size_t locNumParticlesNeedToBeConstrained = 0;
		DKinFitParticle* locConstrainableParticle = NULL;
		deque<DKinFitParticle*> locTempParticles = locKinFitConstraint_P4->dInitialParticles;
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
		locTempParticles = locKinFitConstraint_P4->dFinalParticles;
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
	for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dInitialParticles.size(); ++loc_j)
	{
		locKinFitParticle = locKinFitConstraint_P4->dInitialParticles[loc_j];
		if(locKinFitParticle == locParticleToConstrain)
		{
			locKinFitConstraint_P4->Set_ConstrainedP4Particle(locParticleToConstrain);
			locParticleToConstrain->Set_DecayingParticleAtProductionVertexFlag(false);
			locConstrainedParticleIsInInitialState = true;
		}
		else
			locMomentum += locKinFitParticle->Get_Momentum();
	}
	for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
	{
		locKinFitParticle = locKinFitConstraint_P4->dFinalParticles[loc_j];
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
			if(locKinFitConstraint_P4->dConstrainedP4Particle != dKinFitParticles[loc_i])
				continue; //either no particle p4 constrained, or is the wrong particle
			locP4ConstrainedFlag = true;
			break;
		}
		if(!locP4ConstrainedFlag)
		{
			cout << "ERROR in DKinFitter: Unable to include decaying particle in 2 vertex/time constraints if it's momentum is unknown (not constrained)!!  Exiting." << endl;
			return false; //unable to include decaying particle in vertex/time constraint if it's momentum is unknown (not constrained)
		}

		//set vertex constraint flag of decaying particle
		TVector3 locMomentum = dKinFitParticles[loc_i]->Get_Momentum();
		unsigned short int locVertexConstraintFlag;
		if(fabs(locMomentum.Pz()) > fabs(locMomentum.Px()))
			locVertexConstraintFlag = (fabs(locMomentum.Pz()) > fabs(locMomentum.Py())) ? 1 : 2;
		else
			locVertexConstraintFlag = (fabs(locMomentum.Px()) > fabs(locMomentum.Py())) ? 3 : 2;
		dKinFitParticles[loc_i]->Set_VertexConstraintFlag(locVertexConstraintFlag);

		//below flag true if decaying particle is an initial particle in the p4 constraint it's set in (defined at)
		bool locInitialStateFlag_P4 = !(dKinFitParticles[loc_i]->Get_DecayingParticleAtProductionVertexFlag());

		//loop over the vertex constraints that the decaying particle is involved in
		bool locKnownVertexNotDefinedFlag = true;
		deque<DKinFitConstraint_VertexBase*>::iterator locCommonVertexConstraintIterator = locVertexAndOrTimeConstraints.end();
		TVector3 locDefinedPosition; //the position of where the decaying particle is defined
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
					//the decaying particle is in the same side of the reaction in both constraints: it is defined here
					dKinFitParticles[loc_i]->Set_Position(locKinFitConstraint_Vertex->Get_CommonVertex());
					locKinFitConstraint_Vertex->Add_NoConstrainParticle(dKinFitParticles[loc_i]);
					if(dDebugLevel > 10)
						cout << "q, mass of decaying particle vertex defined by constraint: " << dKinFitParticles[loc_i]->Get_Charge() << ", " << dKinFitParticles[loc_i]->Get_Mass() << endl;
					locKnownVertexNotDefinedFlag = false;
				}
				else
				{
					//the decaying particle is constrained here
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
					//the decaying particle is in the same side of the reaction in both constraints: it is defined here
					dKinFitParticles[loc_i]->Set_Position(locKinFitConstraint_Spacetime->Get_CommonVertex());
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
			if(locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles.size() < 2)
			{
				cout << "ERROR in DKinFitter: not enough particles in spacetime constraint!!  Exiting." << endl;
				return false; //decaying particle constrained to a vertex, but no point on its trajectory is defined (by another vertex constraint)
			}
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			if(locKinFitConstraint_Vertex->dConstrainVertexParticles.size() < 2)
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
	unsigned int locNumP4ConstraintCounter = 0;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			dNumF += (locNumP4ConstraintCounter == 0) ? 4 : 1; //p4/mass
			++locNumP4ConstraintCounter;
			deque<DKinFitParticle*> locInitialParticles = locKinFitConstraint_P4->dInitialParticles;
			deque<DKinFitParticle*> locFinalParticles = locKinFitConstraint_P4->dFinalParticles;
			for(size_t loc_j = 0; loc_j < locInitialParticles.size(); ++loc_j)
			{
				DKinFitParticleType locKinFitParticleType = locInitialParticles[loc_j]->Get_KinFitParticleType();
				if((locKinFitParticleType == d_DecayingParticle) && (locInitialParticles[loc_j]->Get_NumP4Constraints() == 1))
				{
					dNumXi += 3; //p3 //decaying particle included in only one p4 constraint
					break;
				}
			}
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				DKinFitParticleType locKinFitParticleType = locFinalParticles[loc_j]->Get_KinFitParticleType();
				if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locFinalParticles[loc_j]->Get_NumP4Constraints() == 1)))
				{
					dNumXi += 3; //p3 //missing particle or decaying particle included in only one p4 constraint
					break;
				}
			}
			continue;
		}
		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			dNumXi += 3; //v3
			dNumF += 2*locKinFitConstraint_Vertex->dConstrainVertexParticles.size();
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->dConstrainVertexParticles.size(); ++loc_j)
					cout << locKinFitConstraint_Vertex->dConstrainVertexParticles[loc_j]->Get_Charge() << ", " << locKinFitConstraint_Vertex->dConstrainVertexParticles[loc_j]->Get_Mass() << "; ";
				cout << endl;
			}
			continue;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			dNumXi += 4; //v3, t
			deque<DKinFitParticle*> locConstrainSpacetimeParticles = locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles;
			for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
				dNumF += 3;
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of spacetime vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
					cout << locConstrainSpacetimeParticles[loc_j]->Get_Charge() << ", " << locConstrainSpacetimeParticles[loc_j]->Get_Mass() << ";";
				cout << endl;
			}
			if(Get_IsBFieldNearBeamline())
			{
				size_t locNumChargedConstraintParticles = 0;
				size_t locNumDecayingChargedConstraintParticles = 0;
				for(size_t loc_j = 0; loc_j < locConstrainSpacetimeParticles.size(); ++loc_j)
				{
					if(locConstrainSpacetimeParticles[loc_j]->Get_Charge() == 0)
						continue;
					++locNumChargedConstraintParticles;
					if(locConstrainSpacetimeParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle)
						++locNumDecayingChargedConstraintParticles;
				}
				dNumXi += locNumChargedConstraintParticles; //path length (l) for each accelerating particle
				dNumF += locNumChargedConstraintParticles; //extra constraint due to extra unknown (path length (l) for each accelerating particle)
			}
			dNumF += locKinFitConstraint_Spacetime->dOnlyConstrainTimeParticles.size(); //for each neutral shower
			if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
				++dNumF;
		}
	}

	if(dDebugLevel > 10)
		cout << "DKinFitter: Num measurables, unknowns, constraints = " << dNumEta << ", " << dNumXi << ", " << dNumF << endl;
}

void DKinFitter::Resize_Matrices(void)
{
	if(dF.GetNrows() != static_cast<int>(dNumF))
	{
		dF.ResizeTo(dNumF, 1);
		dS.ResizeTo(dNumF, dNumF);
		dS_Inverse.ResizeTo(dNumF, dNumF);
		dLambda.ResizeTo(dNumF, 1);
		dLambda_T.ResizeTo(1, dNumF);
		dF_dEta.ResizeTo(dNumF, dNumEta);
		dF_dEta_T.ResizeTo(dNumEta, dNumF);
		dF_dXi.ResizeTo(dNumF, dNumXi);
		dF_dXi_T.ResizeTo(dNumXi, dNumF);
	}
	else
	{
		if(dF_dEta.GetNcols() != static_cast<int>(dNumEta))
		{
			dF_dEta.ResizeTo(dNumF, dNumEta);
			dF_dEta_T.ResizeTo(dNumEta, dNumF);
		}
		if(dF_dXi.GetNcols() != static_cast<int>(dNumXi))
		{
			dF_dXi.ResizeTo(dNumF, dNumXi);
			dF_dXi_T.ResizeTo(dNumXi, dNumF);
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
		dU.ResizeTo(dNumXi, dNumXi);
		dU_Inverse.ResizeTo(dNumXi, dNumXi);
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
	dF_dEta.Zero();
	dF_dEta_T.Zero();
	dF_dXi.Zero();
	dF_dXi_T.Zero();

	dS.Zero();
	dS_Inverse.Zero();
	dU.Zero();
	dU_Inverse.Zero();

	dVXi->Zero();
	dVEta->Zero();
}

void DKinFitter::Fill_InputMatrices(void)
{
	//fill dY, dEta, dVY, dXi

	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	int locCharge, locParamIndex, locConstraintIndex_Eta, locConstraintIndex_Xi;
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

	//SETUP dXi (with initial guesses) and constraint equation indices
	locParamIndex = 0;
	locConstraintIndex_Eta = 0;
	locConstraintIndex_Xi = 0;
	DKinFitParticle* locConstrainedKinFitParticle;
	deque<DKinFitParticle*> locNoConstrainParticles;
	deque<DKinFitParticle*> locConstrainParticles;
	bool locRFTimeConstrainedFlag = false;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			locKinFitConstraint_P4->Set_FIndex(locConstraintIndex_Eta);
			deque<DKinFitParticle*> locInitialParticles = locKinFitConstraint_P4->dInitialParticles;
			deque<DKinFitParticle*> locFinalParticles = locKinFitConstraint_P4->dFinalParticles;

			for(size_t loc_j = 0; loc_j < locInitialParticles.size(); ++loc_j)
			{
				DKinFitParticleType locKinFitParticleType = locInitialParticles[loc_j]->Get_KinFitParticleType();
				if(locKinFitParticleType == d_BeamParticle)
				{
					locConstraintIndex_Eta += 4; //p4
					break;
				}
				else if((locKinFitParticleType == d_DecayingParticle) && (locInitialParticles[loc_j]->Get_NumP4Constraints() == 1))
				{
					locConstraintIndex_Eta += 4; //p4
					break;
				}
				else
				{
					locConstraintIndex_Eta += 1; //mass
					break;
				}
			}

			locConstrainedKinFitParticle = locKinFitConstraint_P4->dConstrainedP4Particle;
			if(locConstrainedKinFitParticle == NULL)
				continue; //no missing particles or no unknown p3 (p3 is derivable from known p3's (p4 constraint already applied elsewhere))
			locKinFitParticleType = locConstrainedKinFitParticle->Get_KinFitParticleType();
			if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locConstrainedKinFitParticle->Get_NumP4Constraints() == 1)))
			{
				//set initial p3 guess
				TVector3 locMomentum = locConstrainedKinFitParticle->Get_Momentum();
				dXi(locParamIndex, 0) = locMomentum.Px();
				dXi(locParamIndex + 1, 0) = locMomentum.Py();
				dXi(locParamIndex + 2, 0) = locMomentum.Pz();
				locConstrainedKinFitParticle->Set_PxParamIndex(locParamIndex);
				locParamIndex += 3;
			}
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
			locNoConstrainParticles = locKinFitConstraint_Vertex->dNoConstrainParticles;
			for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
			{
				if((locNoConstrainParticles[loc_j]->Get_KinFitParticleType() == d_MissingParticle) || (locNoConstrainParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle))
					locNoConstrainParticles[loc_j]->Set_VxParamIndex(locParamIndex); //not included in fit, but particle vertex is defined by the fit result
			}
			locParamIndex += 3;
			locConstrainParticles = locKinFitConstraint_Vertex->dConstrainVertexParticles;
			for(size_t loc_j = 0; loc_j < locConstrainParticles.size(); ++loc_j)
			{
				locKinFitConstraint_Vertex->Set_FIndex(locConstrainParticles[loc_j], locConstraintIndex_Eta);
				locConstraintIndex_Eta += 2;
			}
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

			locNoConstrainParticles = locKinFitConstraint_Spacetime->dNoConstrainParticles;
			for(size_t loc_j = 0; loc_j < locNoConstrainParticles.size(); ++loc_j)
			{
				if((locNoConstrainParticles[loc_j]->Get_KinFitParticleType() != d_MissingParticle) && (locNoConstrainParticles[loc_j]->Get_KinFitParticleType() != d_DecayingParticle))
					continue;
				locNoConstrainParticles[loc_j]->Set_VxParamIndex(locParamIndex); //not included in fit, but particle vertex is defined by the fit result
				locNoConstrainParticles[loc_j]->Set_TParamIndex(locParamIndex + 3); //not included in fit, but particle time is defined by the fit result
			}
			locParamIndex += 4;

			locConstrainParticles = locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles;
			for(size_t loc_j = 0; loc_j < locConstrainParticles.size(); ++loc_j)
			{
				if(locConstrainParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle)
				{
					locKinFitConstraint_Spacetime->Set_FIndex(locConstrainParticles[loc_j], locConstraintIndex_Xi);
					locConstraintIndex_Xi += 3;
				}
				else
				{
					locKinFitConstraint_Spacetime->Set_FIndex(locConstrainParticles[loc_j], locConstraintIndex_Eta);
					locConstraintIndex_Eta += 3;
				}
			}

			if(Get_IsBFieldNearBeamline())
			{
				for(size_t loc_j = 0; loc_j < locConstrainParticles.size(); ++loc_j)
				{
					if(locConstrainParticles[loc_j]->Get_Charge() == 0)
						continue;
					locConstrainParticles[loc_j]->Set_LParamIndex(locParamIndex);
					dXi(locParamIndex, 0) = locConstrainParticles[loc_j]->Get_PathLength();
					++locParamIndex;
					if(locConstrainParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle)
						++locConstraintIndex_Xi;
					else
						++locConstraintIndex_Eta;
				}
			}
			locConstrainParticles = locKinFitConstraint_Spacetime->dOnlyConstrainTimeParticles;
			for(size_t loc_j = 0; loc_j < locConstrainParticles.size(); ++loc_j) //neutral showers
			{
				locKinFitConstraint_Spacetime->Set_FIndex(locConstrainParticles[loc_j], locConstraintIndex_Eta);
				++locConstraintIndex_Eta; 
			}
			if((dRFMatchedBeamParticle != NULL) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
			{
				locRFTimeConstrainedFlag = true;
				locKinFitConstraint_Spacetime->Set_FIndex(NULL, locConstraintIndex_Eta);
				++locConstraintIndex_Eta;
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
	// translate data from Eta to particles
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

			deque<DKinFitParticle*> locConstrainSpacetimeParticles = locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles;
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

	// calc non-unknown decaying particle momentum (derived from other particles)
		//must do last because assumes all other particle p3's are updated
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 == NULL)
			continue;
		deque<DKinFitParticle*> locInitialParticles = locKinFitConstraint_P4->dInitialParticles;
		DKinFitParticle* locKinFitParticle = locInitialParticles[0];
		if(locKinFitParticle->Get_PxParamIndex() >= 0)
			continue; //already updated
		//locKinFitParticle is now definitely an enclosed decaying particle
		if(locKinFitParticle->Get_DecayingParticleAtProductionVertexFlag())
			locKinFitParticle->Set_Momentum(Calc_DecayingP4_FromDecayProducts(locKinFitConstraint_P4, false).Vect()); //returns p4 at production vertex
		else
			locKinFitParticle->Set_Momentum(Calc_DecayingP4_FromDecayProducts(locKinFitConstraint_P4, true).Vect()); //returns p4 at decay vertex
	}

	if(dRFTimeParamIndex > 0)
		dRFTime = dEta(dRFTimeParamIndex, 0);
}

TLorentzVector DKinFitter::Calc_DecayingP4_FromDecayProducts(DKinFitConstraint_P4* locP4Constraint, bool locDecayMomentumFlag) const
{
	//locDecayMomentumFlag = true to return momentum at the decay vertex, false at the production vertex
	TLorentzVector locP4;

	if(!locDecayMomentumFlag)
	{
		//propagate decaying particle momentum from the decay vertex to the production vertex
		deque<DKinFitParticle*> locInitialParticles = locP4Constraint->dInitialParticles;
		DKinFitParticle* locKinFitParticle = locInitialParticles[0]; //the decaying particle

		size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
		bool locEnoughVertexFitsFlag = (locNumVertexFits == 2);
		int locCharge = locKinFitParticle->Get_Charge();
		bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();

		if(locEnoughVertexFitsFlag && locChargedBFieldFlag)
		{
			TVector3 locPosition = locKinFitParticle->Get_Position();
			TVector3 locDeltaX = locKinFitParticle->Get_CommonVertex() - locPosition;
			TVector3 locBField = Get_BField(locPosition);
			TVector3 locH = locBField.Unit();
			double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
			locP4.SetVect(locP4.Vect() - locDeltaX.Cross(locA*locH));
		}
	}

	deque<DKinFitParticle*> locFinalParticles = locP4Constraint->dFinalParticles;
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locFinalParticles[loc_i];

		size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		bool locEnoughVertexFitsFlag = (locNumVertexFits > 0) && ((locNumVertexFits == 2) || (locKinFitParticleType != d_DecayingParticle));
		int locCharge = locKinFitParticle->Get_Charge();
		bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();

		if(locEnoughVertexFitsFlag && locChargedBFieldFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle))
		{
			TVector3 locPosition = locKinFitParticle->Get_Position();
			TVector3 locDeltaX = locKinFitParticle->Get_CommonVertex() - locPosition;
			TVector3 locBField = Get_BField(locPosition);
			TVector3 locH = locBField.Unit();
			double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
			locP4.SetVect(locP4.Vect() - locDeltaX.Cross(locA*locH));
		}

		if((locFinalParticles[loc_i]->Get_KinFitParticleType() == d_DecayingParticle) && (locFinalParticles[loc_i]->Get_PxParamIndex() == -1))
		{
			//decaying particle whose momentum must also be derived
			deque<DKinFitConstraint_P4*> locP4Constraints = locFinalParticles[loc_i]->Get_P4Constraints();
			for(size_t loc_j = 0; loc_j < locP4Constraints.size(); ++loc_j)
			{
				if(locP4Constraints[loc_j] == locP4Constraint)
					continue;
				locP4 += Calc_DecayingP4_FromDecayProducts(locP4Constraints[loc_j], true);
				break;
			}
		}
		else
			locP4 += locFinalParticles[loc_i]->Get_P4();
	}
	return locP4;
}

void DKinFitter::Calc_dF(void)
{
	dF.Zero();
	dF_dXi.Zero();
	dF_dEta.Zero();
	size_t locFIndex = 0;
	DKinFitParticle* locKinFitParticle;
	bool locIsDecayingFlag = false;
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			if(dDebugLevel > 10)
				cout << "DKinFitter: F index = " << locKinFitConstraint_P4->Get_FIndex() << endl;

			locKinFitParticle = (locKinFitConstraint_P4->dInitialParticles)[0];
			DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
			int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();

			if((locKinFitParticleType == d_BeamParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex >= 0)))
			{
				//initial particle is beam particle or open-ended decaying particle: p4 constraint instead of mass constraint
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dInitialParticles.size(); ++loc_j)
				{
					locKinFitParticle = (locKinFitConstraint_P4->dInitialParticles)[loc_j];
					Calc_dF_P4(locKinFitConstraint_P4, locKinFitParticle, true, false, true, NULL);
				}
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
				{
					locKinFitParticle = (locKinFitConstraint_P4->dFinalParticles)[loc_j];
					Calc_dF_P4(locKinFitConstraint_P4, locKinFitParticle, false, false, false, NULL);
				}
			}
			else
			{
				//mass constraint
				Calc_dF_Mass(locKinFitConstraint_P4);
				TLorentzVector locP4FromDecayProducts = Calc_DecayingP4_FromDecayProducts(locKinFitConstraint_P4, true);
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
				{
					locKinFitParticle = (locKinFitConstraint_P4->dFinalParticles)[loc_j];
					Calc_dF_Mass_Derivs(locKinFitConstraint_P4, locP4FromDecayProducts, locKinFitParticle);
				}
			}
			continue;
		}

		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->dConstrainVertexParticles.size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Vertex->dConstrainVertexParticles)[loc_j];
				DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

				locIsDecayingFlag = (locKinFitParticleType == d_DecayingParticle);
				locFIndex = locKinFitConstraint_Vertex->Get_FIndex(locKinFitParticle);
				if(dDebugLevel > 10)
					cout << "DKinFitter: F index, locIsDecayingFlag = " << locFIndex << ", " << locIsDecayingFlag << endl;

				bool locInitialStateFlag = false; //unless set otherwise
				if(locKinFitParticleType == d_BeamParticle)
					locInitialStateFlag = true;
				//for locDecayingParticles: bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
				deque<pair<DKinFitParticle*, bool> > locDecayingParticles = locKinFitConstraint_Vertex->dDecayingParticles;
				for(size_t loc_k = 0; loc_k < locDecayingParticles.size(); ++loc_k)
				{
					if(locDecayingParticles[loc_k].first != locKinFitParticle)
						continue;
					locInitialStateFlag = !locDecayingParticles[loc_k].second;
					break;
				}

				Calc_dF_Vertex(locFIndex, locKinFitParticle, NULL, locInitialStateFlag, locInitialStateFlag);
			}
		}

		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles.size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Spacetime->dConstrainSpacetimeParticles)[loc_j];
				locIsDecayingFlag = (locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle);
				locFIndex = locKinFitConstraint_Vertex->Get_FIndex(locKinFitParticle);
				if(dDebugLevel > 10)
					cout << "DKinFitter: F index, locIsDecayingFlag = " << locFIndex << ", " << locIsDecayingFlag << endl;
				Calc_dF_Time(locFIndex, locKinFitParticle, false);
				if((locKinFitParticle == dRFMatchedBeamParticle) && locKinFitConstraint_Spacetime->Get_UseRFTimeFlag())
				{
					locFIndex = locKinFitConstraint_Vertex->Get_FIndex(NULL);
					if(dDebugLevel > 10)
						cout << "DKinFitter: F index, locIsDecayingFlag = " << locFIndex << ", " << locIsDecayingFlag << endl;
					Calc_dF_Time(locFIndex, dRFMatchedBeamParticle, true);
				}
			}
		}
	}
	dF_dEta_T.Transpose(dF_dEta);
	dF_dXi_T.Transpose(dF_dXi);
}

void DKinFitter::Calc_dF_P4(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, bool locDerivsOnlyFlag, bool locOriginalInitialStateFlag, DKinFitConstraint_P4* locKinFitSubConstraint_P4)
{
	size_t locFIndex = locKinFitConstraint_P4->Get_FIndex();
	const DKinFitParticle* locConstrainedParticle = locKinFitConstraint_P4->dConstrainedP4Particle;

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

	size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
	bool locEnoughVertexFitsFlag = (locNumVertexFits > 0) && ((locNumVertexFits == 2) || (locKinFitParticleType != d_DecayingParticle));
	bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();
	double locSignMultiplier = locInitialStateFlag ? 1.0 : -1.0;

	if(!locDerivsOnlyFlag)
	{
		dF(locFIndex, 0) += locSignMultiplier*locP4.E();
		dF(locFIndex + 1, 0) += locSignMultiplier*locP4.Px();
		dF(locFIndex + 2, 0) += locSignMultiplier*locP4.Py();
		dF(locFIndex + 3, 0) += locSignMultiplier*locP4.Pz();

		if(locEnoughVertexFitsFlag && locChargedBFieldFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle) && (locKinFitParticle != locConstrainedParticle))
		{
			TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
			//fitting vertex of charged track in magnetic field that is not the constrained (decaying or missing) particle: momentum changes as function of vertex
			dF(locFIndex + 1, 0) -= locSignMultiplier*locA*locDeltaXCrossH.X();
			dF(locFIndex + 2, 0) -= locSignMultiplier*locA*locDeltaXCrossH.Y();
			dF(locFIndex + 3, 0) -= locSignMultiplier*locA*locDeltaXCrossH.Z();
		}
	}

	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 1; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		return; //target params are fixed: no partial derivatives
	}
	else if(locChargedBFieldFlag && locEnoughVertexFitsFlag && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle) & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 2; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dF_dEta(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dF_dEta(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;

		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locA*locH.Z();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locSignMultiplier*locA*locH.Y();

		dF_dEta(locFIndex + 2, locVxParamIndex) = -1.0*locSignMultiplier*locA*locH.Z();
		dF_dEta(locFIndex + 2, locVxParamIndex + 2) = locSignMultiplier*locA*locH.X();

		dF_dEta(locFIndex + 3, locVxParamIndex) = locSignMultiplier*locA*locH.Y();
		dF_dEta(locFIndex + 3, locVxParamIndex + 1) = -1.0*locSignMultiplier*locA*locH.X();

		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);

		dF_dXi(locFIndex + 2, locCommonVxParamIndex) -= dF_dEta(locFIndex + 2, locVxParamIndex);
		dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 2, locVxParamIndex + 2);

		dF_dXi(locFIndex + 3, locCommonVxParamIndex) -= dF_dEta(locFIndex + 3, locVxParamIndex);
		dF_dXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 3, locVxParamIndex + 1);
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 3; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locEParamIndex) = locSignMultiplier;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dF_dEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Px();
		dF_dEta(locFIndex + 2, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Py();
		dF_dEta(locFIndex + 3, locEParamIndex) = locSignMultiplier*locEOverPSq*locP4.Pz();

		TVector3 locDeltaXOverMagDeltaXSq = locDeltaX*(1.0/locDeltaX.Mag2());

		dF_dEta(locFIndex + 1, locVxParamIndex) = locSignMultiplier*locP4.Px()*(locDeltaXOverMagDeltaXSq.X() - 1.0/locDeltaX.X());
		dF_dEta(locFIndex + 2, locVxParamIndex + 1) = locSignMultiplier*locP4.Py()*(locDeltaXOverMagDeltaXSq.Y() - 1.0/locDeltaX.Y());
		dF_dEta(locFIndex + 3, locVxParamIndex + 2) = locSignMultiplier*locP4.Pz()*(locDeltaXOverMagDeltaXSq.Z() - 1.0/locDeltaX.Z());

		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Y();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Z();

		dF_dEta(locFIndex + 2, locVxParamIndex) = locSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.X();
		dF_dEta(locFIndex + 2, locVxParamIndex + 2) = locSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.Z();

		dF_dEta(locFIndex + 3, locVxParamIndex) = locSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.X();
		dF_dEta(locFIndex + 3, locVxParamIndex + 1) = locSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.Y();

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 2, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 2, locVxParamIndex + 1);
		dF_dXi(locFIndex + 3, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 3, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);

		dF_dXi(locFIndex + 2, locCommonVxParamIndex) -= dF_dEta(locFIndex + 2, locVxParamIndex);
		dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 2, locVxParamIndex + 2);

		dF_dXi(locFIndex + 3, locCommonVxParamIndex) -= dF_dEta(locFIndex + 3, locVxParamIndex);
		dF_dXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 3, locVxParamIndex + 1);
	}
	else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex >= 0)))
	{
		//missing or open-ended-decaying particle: p3 is unknown (not derivable) //must be the constrained particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 4; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dXi(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dF_dXi(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dF_dXi(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dF_dXi(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dF_dXi(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dF_dXi(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;
	}
	else if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 5; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticle->Get_NumVertexFits() == 2))
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_P4() Section 5a" << endl;
			//vertex factors
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locSignMultiplier*locA*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += -1.0*locSignMultiplier*locA*locH.Y();

			dF_dXi(locFIndex + 2, locVxParamIndex) += -1.0*locSignMultiplier*locA*locH.Z();
			dF_dXi(locFIndex + 2, locVxParamIndex + 2) += locSignMultiplier*locA*locH.X();

			dF_dXi(locFIndex + 3, locVxParamIndex) += locSignMultiplier*locA*locH.Y();
			dF_dXi(locFIndex + 3, locVxParamIndex + 1) += -1.0*locSignMultiplier*locA*locH.X();

			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);

			dF_dXi(locFIndex + 2, locCommonVxParamIndex) -= dF_dXi(locFIndex + 2, locVxParamIndex);
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 2, locVxParamIndex + 2);

			dF_dXi(locFIndex + 3, locCommonVxParamIndex) -= dF_dXi(locFIndex + 3, locVxParamIndex);
			dF_dXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 3, locVxParamIndex + 1);
		}

		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->Get_P4Constraints();
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		{
			DKinFitConstraint_P4* locNewKinFitSubConstraint_P4 = locP4Constraints[loc_i];
			if(locNewKinFitSubConstraint_P4 == locKinFitConstraint_P4)
				continue;
			if(locNewKinFitSubConstraint_P4 == locKinFitSubConstraint_P4)
				continue;

			//if this constraint contains the decaying particle, replace it with the other particles it's momentum is derived from
			deque<DKinFitParticle*> locInitialParticles = locNewKinFitSubConstraint_P4->dInitialParticles;
			deque<DKinFitParticle*> locFinalParticles = locNewKinFitSubConstraint_P4->dFinalParticles;
			for(size_t loc_j = 0; loc_j < locInitialParticles.size(); ++loc_j)
			{
				if(locInitialParticles[loc_j] == locKinFitParticle)
					continue;
//detected: +/- if init/final
//p-replacement: factor of +/- if same/different state
//so if ORIG decaying was init/final (+/-), and detected is init (+), then should call as if it was (+/+): init
//so if ORIG decaying was init/final (+/-), and detected is final (-), then should call as if it was (-/-): final
				if(dDebugLevel > 30)
					cout << "decaying, partially replace with init-state q, mass = " << locInitialParticles[loc_j]->Get_Charge() << ", " << locInitialParticles[loc_j]->Get_Mass() << endl;
				Calc_dF_P4(locKinFitConstraint_P4, locInitialParticles[loc_j], true, true, locOriginalInitialStateFlag, locNewKinFitSubConstraint_P4); //else !locInitialStateFlag
			}
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				if(locFinalParticles[loc_j] == locKinFitParticle)
					continue;
				if(dDebugLevel > 30)
					cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
				Calc_dF_P4(locKinFitConstraint_P4, locFinalParticles[loc_j], false, true, locOriginalInitialStateFlag, locNewKinFitSubConstraint_P4);
			}
		}
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 6; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = locSignMultiplier*locP4.Px()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locP4.Py()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locP4.Pz()/locP4.E();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier;
		dF_dEta(locFIndex + 2, locPxParamIndex + 1) = locSignMultiplier;
		dF_dEta(locFIndex + 3, locPxParamIndex + 2) = locSignMultiplier;
	}
}

void DKinFitter::Calc_dF_Mass(DKinFitConstraint_P4* locKinFitConstraint_P4)
{
	size_t locFIndex = locKinFitConstraint_P4->Get_FIndex();

	deque<DKinFitParticle*> locFinalParticles = locKinFitConstraint_P4->dFinalParticles;
	TLorentzVector locP4(0.0, 0.0, 0.0, 0.0);
	double locTargetedMass = (locKinFitConstraint_P4->dInitialParticles)[0]->Get_Mass();

	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locFinalParticles[loc_i];
		locP4 += locKinFitParticle->Get_P4();

		size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		bool locEnoughVertexFitsFlag = (locNumVertexFits > 0) && ((locNumVertexFits == 2) || (locKinFitParticleType != d_DecayingParticle));
		int locCharge = locKinFitParticle->Get_Charge();
		bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();

		if(locEnoughVertexFitsFlag && locChargedBFieldFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle))
		{
			TVector3 locPosition = locKinFitParticle->Get_Position();
			TVector3 locDeltaX = locKinFitParticle->Get_CommonVertex() - locPosition;
			TVector3 locBField = Get_BField(locPosition);
			TVector3 locH = locBField.Unit();
			double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
			locP4.SetVect(locP4.Vect() - locDeltaX.Cross(locA*locH));
		}

	}
	dF(locFIndex, 0) += locP4.M2() - locTargetedMass*locTargetedMass;
}

void DKinFitter::Calc_dF_Mass_Derivs(DKinFitConstraint_P4* locKinFitConstraint_P4, const TLorentzVector& locDecayingP4_FromDecayProducts, const DKinFitParticle* locKinFitParticle)
{
	size_t locFIndex = locKinFitConstraint_P4->Get_FIndex();

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

	size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
	bool locEnoughVertexFitsFlag = (locNumVertexFits > 0) && ((locNumVertexFits == 2) || (locKinFitParticleType != d_DecayingParticle));
	bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	TVector3 locPXCrossHi = locDecayingP4_FromDecayProducts.Vect().Cross(locH);

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 1; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		return; //target params are fixed: no partial derivatives
	}
	else if(locChargedBFieldFlag && locEnoughVertexFitsFlag && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle) & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 2; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*(locP4.Px()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*(locP4.Py()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*(locP4.Pz()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Pz());

		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locA*locPXCrossHi.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locA*locPXCrossHi.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locA*locPXCrossHi.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 3; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dF_dEta(locFIndex, locEParamIndex) = 2.0*(locDecayingP4_FromDecayProducts.E() - locEOverPSq*locDecayingP4_FromDecayProducts.Vect().Dot(locP4.Vect()));

		double locDeltaXDotPXOverMagDeltaXSq = locDeltaX.Dot(locDecayingP4_FromDecayProducts.Vect())/(locDeltaX.Mag2());
		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locP4.Px()*(locDecayingP4_FromDecayProducts.Px()/locDeltaX.X() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locP4.Py()*(locDecayingP4_FromDecayProducts.Py()/locDeltaX.Y() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locP4.Pz()*(locDecayingP4_FromDecayProducts.Pz()/locDeltaX.Z() - locDeltaXDotPXOverMagDeltaXSq);

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex >= 0)))
	{
		//missing or open-ended-decaying particle: p3 is unknown (not derivable) //must be the constrained particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 4; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		dF_dXi(locFIndex, locPxParamIndex) = 2.0*(locP4.Px()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Px());
		dF_dXi(locFIndex, locPxParamIndex + 1) = 2.0*(locP4.Py()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Py());
		dF_dXi(locFIndex, locPxParamIndex + 2) = 2.0*(locP4.Pz()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Pz());
	}
	else if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 5; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticle->Get_NumVertexFits() == 2))
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 5a" << endl;

			dF_dXi(locFIndex, locVxParamIndex) += 2.0*locA*locPXCrossHi.X();
			dF_dXi(locFIndex, locVxParamIndex + 1) += 2.0*locA*locPXCrossHi.Y();
			dF_dXi(locFIndex, locVxParamIndex + 2) += 2.0*locA*locPXCrossHi.Z();

			dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);
		}
		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->Get_P4Constraints();
		DKinFitConstraint_P4* locKinFitSubConstraint_P4 = NULL;
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		{
			if(locP4Constraints[loc_i] == locKinFitConstraint_P4)
				continue;
			locKinFitSubConstraint_P4 = locP4Constraints[loc_i];
			break;
		}
		if(locKinFitSubConstraint_P4 == NULL)
			return; //this shouldn't be possible...
		//ok, now replace the contribution of the decay particle with that of the other particles it's momentum is derived from
		deque<DKinFitParticle*> locFinalParticles = locKinFitSubConstraint_P4->dFinalParticles;
		for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
		{
			if(locFinalParticles[loc_j] == locKinFitParticle)
				continue;
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
			Calc_dF_Mass_Derivs(locKinFitConstraint_P4, locDecayingP4_FromDecayProducts, locFinalParticles[loc_j]);
		}
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Mass_Derivs() Section 6; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*(locP4.Px()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*(locP4.Py()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*(locP4.Pz()*locDecayingP4_FromDecayProducts.E()/locP4.E() - locDecayingP4_FromDecayProducts.Pz());
	}
}

void DKinFitter::Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag)
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 1; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		return; //no partial derivatives
	}

	if(locKinFitParticle_DecayingSource == NULL)
		Calc_dF_Vertex_NotDecaying(locFIndex, locKinFitParticle);
	else
	{
		int locCharge = locKinFitParticle_DecayingSource->Get_Charge();
		if((locCharge != 0) && Get_IsBFieldNearBeamline())
			Calc_dF_Vertex_Decaying_Accel(locFIndex, locKinFitParticle, locKinFitParticle_DecayingSource, locInitialStateFlag, locOriginalInitialStateFlag);
		else
			Calc_dF_Vertex_Decaying_NonAccel(locFIndex, locKinFitParticle, locKinFitParticle_DecayingSource, locInitialStateFlag, locOriginalInitialStateFlag);
	}

	if(locKinFitParticleType == d_DecayingParticle)
	{
		//replace the momentum with it's component tracks //call this func for each decay product
		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->Get_P4Constraints();
		DKinFitConstraint_P4* locKinFitConstraint_P4 = NULL;
		// find the p4 constraint where the decaying particle momentum is constrained
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		{
			if(locP4Constraints[loc_i]->dConstrainedP4Particle != locKinFitParticle)
				continue;
			locKinFitConstraint_P4 = locP4Constraints[loc_i];
			break;
		}
		if(locKinFitConstraint_P4 == NULL)
			return; //shouldn't be possible...

		//ok, now replace the contribution of the decay particle with that of the other particles it's momentum is derived from
		deque<DKinFitParticle*> locInitialParticles = locKinFitConstraint_P4->dInitialParticles;
		deque<DKinFitParticle*> locFinalParticles = locKinFitConstraint_P4->dFinalParticles;
		for(size_t loc_j = 0; loc_j < locInitialParticles.size(); ++loc_j)
		{
			if(locInitialParticles[loc_j] == locKinFitParticle)
				continue;
//detected: +/- if init/final
//p-replacement: factor of +/- if same/different state
//so if ORIG decaying was init/final (+/-), and detected is init (+), then should call as if it was (+/+): init
//so if ORIG decaying was init/final (+/-), and detected is final (-), then should call as if it was (-/-): final
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state q, mass = " << locInitialParticles[loc_j]->Get_Charge() << ", " << locInitialParticles[loc_j]->Get_Mass() << endl;
			Calc_dF_Vertex(locFIndex, locInitialParticles[loc_j], locKinFitParticle, true, locOriginalInitialStateFlag);
		}
		for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
		{
			if(locFinalParticles[loc_j] == locKinFitParticle)
				continue;
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
			Calc_dF_Vertex(locFIndex, locFinalParticles[loc_j], locKinFitParticle, false, locOriginalInitialStateFlag);
		}
	}
}

void DKinFitter::Calc_dF_Vertex_NotDecaying(size_t locFIndex, const DKinFitParticle* locKinFitParticle)
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

	if((locCharge != 0) && Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle))) //DONE
	{
		//detected charged particle in b-field (can be beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 2; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locBField = Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		TVector3 locPCrossH = locMomentum.Cross(locH);
		TVector3 locPCrossDeltaX = locMomentum.Cross(locDeltaX);
		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);

		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle, locJ, locQ, locM, locD);

		dF(locFIndex, 0) = locPCrossDeltaX.Dot(locH) - 0.5*locA*(locDeltaX.Mag2() - locDeltaXDotH*locDeltaXDotH);
		dF(locFIndex + 1, 0) = locDeltaXDotH - (locPDotH/locA)*asin(locJ);

		dF_dEta(locFIndex, locPxParamIndex) = locDeltaXCrossH.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaXCrossH.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locDeltaXCrossH.Z();

		dF_dEta(locFIndex, locVxParamIndex) = (locPCrossH + locA*locM).X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = (locPCrossH + locA*locM).Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = (locPCrossH + locA*locM).Z();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locQ.X();
		dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locQ.Y();
		dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locQ.Z();

		dF_dEta(locFIndex + 1, locVxParamIndex) = locD.X();
		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locD.Y();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locD.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle)) //DONE
	{
		//constraining this decaying charged particle in b-field //one-time contributions from decaying particle, does not include the particles it is replaced by (elsewhere)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 3; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle, locJ, locQ, locM, locD);

		TVector3 locBField = Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locPCrossH = locMomentum.Cross(locH);
		TVector3 locPCrossDeltaX = locMomentum.Cross(locDeltaX);
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);

		dF(locFIndex, 0) = locPCrossDeltaX.Dot(locH) - 0.5*locA*(locDeltaX.Mag2() - locDeltaXDotH*locDeltaXDotH);
		dF(locFIndex + 1, 0) = locDeltaXDotH - (locPDotH/locA)*asin(locJ);

		dF_dXi(locFIndex, locVxParamIndex) = locPCrossH.X() + locA*locM.X();
		dF_dXi(locFIndex, locVxParamIndex + 1) = locPCrossH.Y() + locA*locM.Y();
		dF_dXi(locFIndex, locVxParamIndex + 2) = locPCrossH.Z() + locA*locM.Z();

		dF_dXi(locFIndex + 1, locVxParamIndex) = locD.X();
		dF_dXi(locFIndex + 1, locVxParamIndex + 1) = locD.Y();
		dF_dXi(locFIndex + 1, locVxParamIndex + 2) = locD.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dXi(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticleType == d_DecayingParticle) //non-accel decaying particle //DONE
	{
		//constraining this decaying non-accel particle //one-time contributions from decaying particle, does not include the particles it is replaced by (elsewhere)

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 2);
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		}
		else //2 & 3 //px is largest
		{
			dF(locFIndex, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) = locMomentum.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex) = -1.0*dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		}
	}
	else //neutral detected or beam particle, or no magnetic field //DONE
	{
		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 2);
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		}
		else //2 & 3 //px is largest
		{
			dF(locFIndex, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locDeltaX.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locMomentum.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		}
	}
}

void DKinFitter::Calc_dF_Vertex_Decaying_Accel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag)
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
	int locCharge = locKinFitParticle->Get_Charge();

	double locEnergy = locKinFitParticle->Get_P4().E();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locCharge != 0) && Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle))) //DONE
	{
		//detected charged particle in b-field (can be beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 6; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;
		TVector3 locMomentum_DecayingSource = locKinFitParticle_DecayingSource->Get_Momentum();

		TVector3 locBField_DecayingSource = Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();
		TVector3 locBField = Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		TVector3 locDeltaXCrossH_DecayingSource_CrossH = locDeltaXCrossH_DecayingSource.Cross(locH);
		TVector3 locQCrossH = locQ.Cross(locH);

		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		dF_dEta(locFIndex, locPxParamIndex) = locSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dEta(locFIndex, locVxParamIndex) = -1.0*locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Z();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier*locQ.X();
		dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locSignMultiplier*locQ.Y();
		dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locSignMultiplier*locQ.Z();

		dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locSignMultiplier*locA*locQCrossH.X();
		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locSignMultiplier*locA*locQCrossH.Y();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locSignMultiplier*locA*locQCrossH.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag()) //DONE
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 7; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;
		TVector3 locMomentum_DecayingSource = locKinFitParticle_DecayingSource->Get_Momentum();

		TVector3 locBField_DecayingSource = Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		double locXTerm = locDeltaXCrossH_DecayingSource.Dot(locDeltaX)/(locDeltaX.Mag2());
		double locQTerm = locQ.Dot(locDeltaX)/(locDeltaX.Mag2());

		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		dF_dEta(locFIndex, locEParamIndex) = locSignMultiplier*locEnergy*(locDeltaXCrossH_DecayingSource.Dot(locMomentum))/(locMomentum.Mag2());
		dF_dEta(locFIndex, locVxParamIndex) = -1.0*locSignMultiplier*locMomentum.Px()*(locDeltaXCrossH_DecayingSource.X()/(locDeltaX.X()) - locXTerm);
		dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locSignMultiplier*locMomentum.Py()*(locDeltaXCrossH_DecayingSource.Y()/(locDeltaX.Y()) - locXTerm);
		dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locSignMultiplier*locMomentum.Pz()*(locDeltaXCrossH_DecayingSource.Z()/(locDeltaX.Z()) - locXTerm);

		dF_dEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEnergy*(locMomentum.Dot(locQ))/(locMomentum.Mag2());
		dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locSignMultiplier*locMomentum.Px()*(locQ.X()/(locDeltaX.X()) - locQTerm);
		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locSignMultiplier*locMomentum.Py()*(locQ.Y()/(locDeltaX.Y()) - locQTerm);
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locSignMultiplier*locMomentum.Pz()*(locQ.Z()/(locDeltaX.Z()) - locQTerm);

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)) //DONE
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 8; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		//detected/beam, non-accelerating particle
		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;
		TVector3 locMomentum_DecayingSource = locKinFitParticle_DecayingSource->Get_Momentum();

		TVector3 locBField_DecayingSource = Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);

		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		dF_dEta(locFIndex, locPxParamIndex) = locSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier*locQ.X();
		dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locSignMultiplier*locQ.Y();
		dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locSignMultiplier*locQ.Z();
	}
	else if(locKinFitParticleType == d_MissingParticle) //DONE
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 9; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		//missing particle
		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;
		TVector3 locMomentum_DecayingSource = locKinFitParticle_DecayingSource->Get_Momentum();

		TVector3 locBField_DecayingSource = Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);

		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		dF_dXi(locFIndex, locPxParamIndex) = locSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dXi(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dXi(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dXi(locFIndex + 1, locPxParamIndex) = locSignMultiplier*locQ.X();
		dF_dXi(locFIndex + 1, locPxParamIndex + 1) = locSignMultiplier*locQ.Y();
		dF_dXi(locFIndex + 1, locPxParamIndex + 2) = locSignMultiplier*locQ.Z();
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle)) //DONE
	{
		//decaying charged particle in b-field (doesn't include contributions from the particles it is replaced by here)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 10; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;
		TVector3 locMomentum_DecayingSource = locKinFitParticle_DecayingSource->Get_Momentum();

		TVector3 locBField_DecayingSource = Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();
		TVector3 locBField = Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		TVector3 locDeltaXCrossH_DecayingSource_CrossH = locDeltaXCrossH_DecayingSource.Cross(locH);
		TVector3 locQCrossH = locQ.Cross(locH);

		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		dF_dXi(locFIndex, locVxParamIndex) -= locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.X();
		dF_dXi(locFIndex, locVxParamIndex + 1) -= locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Y();
		dF_dXi(locFIndex, locVxParamIndex + 2) -= locSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Z();

		dF_dXi(locFIndex + 1, locVxParamIndex) -= locSignMultiplier*locA*locQCrossH.X();
		dF_dXi(locFIndex + 1, locVxParamIndex + 1) -= locSignMultiplier*locA*locQCrossH.Y();
		dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locSignMultiplier*locA*locQCrossH.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dXi(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);
	}
}

void DKinFitter::Calc_dF_Vertex_Decaying_NonAccel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag)
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
	int locCharge = locKinFitParticle->Get_Charge();

	double locEnergy = locKinFitParticle->Get_P4().E();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locCharge != 0) && Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle))) //DONE
	{
		//detected charged particle in b-field (can be beam particle)
		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;
		TVector3 locBField = Get_BField(locPosition);
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		TVector3 locH = locBField.Unit();

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 1) = locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 1) = locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locSignMultiplier*locDeltaX_DecayingSource.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag()) //DONE
	{
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locPiCrossDeltaXX = locMomentum.Cross(locDeltaX_DecayingSource);
		TVector3 locDeltaXCrossDeltaXX = locDeltaX.Cross(locDeltaX_DecayingSource);
		double locDeltaXiMagSq = locDeltaX.Mag2();
		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.X()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.Y()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq - locMomentum.Y()*locDeltaX_DecayingSource.Z()/locDeltaX.Y();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq + locMomentum.Z()*locDeltaX_DecayingSource.Y()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq + locMomentum.X()*locDeltaX_DecayingSource.Z()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq;
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq - locMomentum.Z()*locDeltaX_DecayingSource.X()/locDeltaX.Z();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.X()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.Z()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq - locMomentum.Y()*locDeltaX_DecayingSource.Z()/locDeltaX.Y();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq + locMomentum.Z()*locDeltaX_DecayingSource.Y()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq - locMomentum.X()*locDeltaX_DecayingSource.Y()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq + locMomentum.Y()*locDeltaX_DecayingSource.X()/locDeltaX.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq;
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.Y()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locSignMultiplier*locEnergy*locPiCrossDeltaXX.Z()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq + locMomentum.X()*locDeltaX_DecayingSource.Z()/locDeltaX.X();
			dF_dEta(locFIndex, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq - locMomentum.Z()*locDeltaX_DecayingSource.X()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq - locMomentum.X()*locDeltaX_DecayingSource.Y()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq + locMomentum.Y()*locDeltaX_DecayingSource.X()/locDeltaX.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq;
		}
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)) //DONE
	{
		//detected/beam, non-accelerating particle
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
	}
	else if(locKinFitParticleType == d_MissingParticle)
	{
		//missing
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex + 1, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
			dF_dXi(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle)) //DONE
	{
		//decaying charged particle in b-field (doesn't include contributions from the particles it is replaced by here)
		double locSignMultiplier = (locInitialStateFlag == locOriginalInitialStateFlag) ? 1.0 : -1.0;
		TVector3 locBField = Get_BField(locPosition);
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		TVector3 locH = locBField.Unit();

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15a; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 1) += locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15b; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 1) += locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());

		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15c; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locVxParamIndex) += locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dXi(locFIndex, locVxParamIndex + 1) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locA*locSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locA*locSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dXi(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);
	}
}

void DKinFitter::Calc_Vertex_Params(const DKinFitParticle* locKinFitParticle, double& locJ, TVector3& locQ, TVector3& locM, TVector3& locD)
{
	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

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
	locJ = locA*locK/locPCrossHMagSq;
	double locC = locPDotH/(locPCrossHMagSq*sqrt(1.0 - locJ*locJ));

	TVector3 locPCrossHCrossH = locPCrossH.Cross(locH);

	locM = locDeltaX - locDeltaXDotH*locH;
	locD = locC*(locMomentum - locPDotH*locH) - locH;
	locQ = -1.0*locH*(asin(locJ)/locA) - locC*(locM + 2.0*(locK/locPCrossHMagSq)*locPCrossHCrossH);
}

void DKinFitter::Calc_dF_Time(size_t locFIndex, const DKinFitParticle* locKinFitParticle, bool locUseRFTimeFlag)
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
			cout << "DKinFitter: Calc_dF_Time() Section 1; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF(locFIndex, 0) = locDeltaT - locPathLength*locEOverC/locMomentum.Mag();

		TVector3 locMomentumTerm = locMomentum*locMass*locMass*locPathLength*(1.0/(29.9792458*locEnergy*locMomentum.Mag()*locMomentum.Mag2()));
		dF_dEta(locFIndex, locPxParamIndex) = locMomentumTerm.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locMomentumTerm.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locMomentumTerm.Z();
		dF_dEta(locFIndex, locTParamIndex) = -1.0;

		dF_dXi(locFIndex, locCommonTParamIndex) = 1.0;
		dF_dXi(locFIndex, locLParamIndex) = -1.0*locEOverC/locMomentum.Mag();
	}
	else if(locUseRFTimeFlag) //beam particle: neutral or no b-field
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Time() Section 2; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		double locPDotDeltaX = locDeltaX.Dot(locMomentum);
		double locPMagSq = locMomentum.Mag2();

		dF(locFIndex, 0) = locDeltaT - locEOverC*locPDotDeltaX/locPMagSq;

		TVector3 locW = locMomentum*locPDotDeltaX*(locEnergy*locEnergy + locMass*locMass)*(1.0/(29.9792458*locEnergy*locPMagSq*locPMagSq)) - locDeltaX*(locEOverC/locPMagSq);

		dF_dEta(locFIndex, locPxParamIndex) = locW.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locW.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locW.Z();

		TVector3 locVertexTerm = locMomentum*(locEOverC/locPMagSq);

		dF_dEta(locFIndex, locVxParamIndex) = locVertexTerm.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = locVertexTerm.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = locVertexTerm.Z();
		dF_dEta(locFIndex, locTParamIndex) = -1.0;

		dF_dXi(locFIndex, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
		dF_dXi(locFIndex, locCommonTParamIndex) = 1.0;
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag()) //neutral shower
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Time() Section 3; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF(locFIndex, 0) = locDeltaT + locShowerEnergy*locDeltaX.Mag()/(29.9792458*locMomentum.Mag());

		dF_dEta(locFIndex, locEParamIndex) = -1.0*locMass*locMass*locDeltaX.Mag()/(29.9792458*locMomentum.Mag()*locMomentum.Mag2());

		TVector3 locPositionTerm = locDeltaX.Unit()*(locEOverC/(locMomentum.Mag()));

		dF_dEta(locFIndex, locVxParamIndex) = -1.0*locPositionTerm.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locPositionTerm.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locPositionTerm.Z();
		dF_dEta(locFIndex, locTParamIndex) = -1.0;

		dF_dXi(locFIndex, locCommonVxParamIndex) = locPositionTerm.X();
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) = locPositionTerm.Y();
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) = locPositionTerm.Z();
		dF_dXi(locFIndex, locCommonTParamIndex) = 1.0;
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
				cout << "DKinFitter: Calc_dF_Time() Section 4; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locPxParamIndex) = locG.X();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = locG.Y();
			dF_dXi(locFIndex + 2, locPxParamIndex + 2) = locG.Z();

			dF_dXi(locFIndex, locPxParamIndex + 1) = locDkx.Y() + locHlocCOverA.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = locDkx.Z() - locHlocCOverA.Y();

			dF_dXi(locFIndex + 1, locPxParamIndex) = locDky.X() - locHlocCOverA.Z();
			dF_dXi(locFIndex + 1, locPxParamIndex + 2) = locDky.Z() + locHlocCOverA.X();

			dF_dXi(locFIndex + 2, locPxParamIndex) = locDkz.X() + locHlocCOverA.Y();
			dF_dXi(locFIndex + 2, locPxParamIndex + 1) = locDkz.Y() - locHlocCOverA.X();

			dF_dXi(locFIndex + 3, locPxParamIndex) = locZ*locMomentum.X();
			dF_dXi(locFIndex + 3, locPxParamIndex + 1) = locZ*locMomentum.Y();
			dF_dXi(locFIndex + 3, locPxParamIndex + 2) = locZ*locMomentum.Z();

			dF_dXi(locFIndex, locVxParamIndex) = -1.0;
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) = -1.0;
			dF_dXi(locFIndex + 2, locVxParamIndex + 2) = -1.0;
			dF_dXi(locFIndex + 3, locTParamIndex) = -1.0;

			dF_dXi(locFIndex, locCommonVxParamIndex) = 1.0;
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = 1.0;
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) = 1.0;
			dF_dXi(locFIndex + 3, locCommonTParamIndex) = 1.0;

			dF_dXi(locFIndex, locLParamIndex) = locR.X();
			dF_dXi(locFIndex + 1, locLParamIndex + 1) = locR.Y();
			dF_dXi(locFIndex + 2, locLParamIndex + 2) = locR.Z();
			dF_dXi(locFIndex + 3, locLParamIndex) = -1.0*locEOverC/locPMag;
		}
		else //detected particle or beam particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 5; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = locG.X();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locG.Y();
			dF_dEta(locFIndex + 2, locPxParamIndex + 2) = locG.Z();

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDkx.Y() + locHlocCOverA.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locDkx.Z() - locHlocCOverA.Y();

			dF_dEta(locFIndex + 1, locPxParamIndex) = locDky.X() - locHlocCOverA.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locDky.Z() + locHlocCOverA.X();

			dF_dEta(locFIndex + 2, locPxParamIndex) = locDkz.X() + locHlocCOverA.Y();
			dF_dEta(locFIndex + 2, locPxParamIndex + 1) = locDkz.Y() - locHlocCOverA.X();

			dF_dEta(locFIndex + 3, locPxParamIndex) = locZ*locMomentum.X();
			dF_dEta(locFIndex + 3, locPxParamIndex + 1) = locZ*locMomentum.Y();
			dF_dEta(locFIndex + 3, locPxParamIndex + 2) = locZ*locMomentum.Z();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0;
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0;
			dF_dEta(locFIndex + 2, locVxParamIndex + 2) = -1.0;
			dF_dEta(locFIndex + 3, locTParamIndex) = -1.0;

			dF_dXi(locFIndex, locCommonVxParamIndex) = 1.0;
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = 1.0;
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) = 1.0;
			dF_dXi(locFIndex + 3, locCommonTParamIndex) = 1.0;

			dF_dXi(locFIndex, locLParamIndex) = locR.X();
			dF_dXi(locFIndex + 1, locLParamIndex + 1) = locR.Y();
			dF_dXi(locFIndex + 2, locLParamIndex + 2) = locR.Z();
			dF_dXi(locFIndex + 3, locLParamIndex) = -1.0*locEOverC/locPMag;
		}
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
				cout << "DKinFitter: Calc_dF_Time() Section 6; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dXi(locFIndex, locPxParamIndex) = locEOverCP.X()*locDeltaXOverPComponents.X() - locDeltaXOverCE.X();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = locEOverCP.Y()*locDeltaXOverPComponents.Y() - locDeltaXOverCE.Y();
			dF_dXi(locFIndex + 2, locPxParamIndex + 2) = locEOverCP.Z()*locDeltaXOverPComponents.Z() - locDeltaXOverCE.Z();

			dF_dXi(locFIndex, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Y();
			dF_dXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Z();

			dF_dXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.X();
			dF_dXi(locFIndex + 1, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.Z();

			dF_dXi(locFIndex + 2, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.X();
			dF_dXi(locFIndex + 2, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.Y();

			dF_dXi(locFIndex, locVxParamIndex) = locEOverCP.X();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) = locEOverCP.Y();
			dF_dXi(locFIndex + 2, locVxParamIndex + 2) = locEOverCP.Z();

			dF_dXi(locFIndex, locTParamIndex) = -1.0;
			dF_dXi(locFIndex + 1, locTParamIndex) = -1.0;
			dF_dXi(locFIndex + 2, locTParamIndex) = -1.0;

			dF_dXi(locFIndex, locCommonVxParamIndex) = -1.0*dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 1);
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) = -1.0*dF_dXi(locFIndex + 2, locVxParamIndex + 2);

			dF_dXi(locFIndex, locCommonTParamIndex) = 1.0;
			dF_dXi(locFIndex + 1, locCommonTParamIndex) = 1.0;
			dF_dXi(locFIndex + 2, locCommonTParamIndex) = 1.0;
		}
		else //detected particle or beam particle
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Time() Section 7; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = locEOverCP.X()*locDeltaXOverPComponents.X() - locDeltaXOverCE.X();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locEOverCP.Y()*locDeltaXOverPComponents.Y() - locDeltaXOverCE.Y();
			dF_dEta(locFIndex + 2, locPxParamIndex + 2) = locEOverCP.Z()*locDeltaXOverPComponents.Z() - locDeltaXOverCE.Z();

			dF_dEta(locFIndex, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Y();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.X()*locPOverCE.Z();

			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.X();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = -1.0*locDeltaXOverPComponents.Y()*locPOverCE.Z();

			dF_dEta(locFIndex + 2, locPxParamIndex) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.X();
			dF_dEta(locFIndex + 2, locPxParamIndex + 1) = -1.0*locDeltaXOverPComponents.Z()*locPOverCE.Y();

			dF_dEta(locFIndex, locVxParamIndex) = locEOverCP.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locEOverCP.Y();
			dF_dEta(locFIndex + 2, locVxParamIndex + 2) = locEOverCP.Z();

			dF_dEta(locFIndex, locTParamIndex) = -1.0;
			dF_dEta(locFIndex + 1, locTParamIndex) = -1.0;
			dF_dEta(locFIndex + 2, locTParamIndex) = -1.0;

			dF_dXi(locFIndex, locCommonVxParamIndex) = -1.0*dF_dEta(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) = -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 1);
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) = -1.0*dF_dEta(locFIndex + 2, locVxParamIndex + 2);

			dF_dXi(locFIndex, locCommonTParamIndex) = 1.0;
			dF_dXi(locFIndex + 1, locCommonTParamIndex) = 1.0;
			dF_dXi(locFIndex + 2, locCommonTParamIndex) = 1.0;
		}
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
			//assume that the correlations remain unchanged: the changes in the parameters and their uncertainties should be small wrst their values
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		locKinFitParticleType = dKinFitParticles[loc_i]->Get_KinFitParticleType();
		if(locKinFitParticleType == d_TargetParticle)
			continue;

		bool locReconstructedParticleFlag = ((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle));
		TMatrixDSym& locKinFitMatrix = locReconstructedParticleFlag ? *dVXi : *dVEta;
		if(locReconstructedParticleFlag)
		{
			TMatrixDSym* locCovarianceMatrix = Get_MatrixDSymResource();
			locCovarianceMatrix->ResizeTo(7, 7);
			locCovarianceMatrix->Zero();
			dKinFitParticles[loc_i]->Set_CovarianceMatrix(locCovarianceMatrix);
		}
		TMatrixDSym& locCovarianceMatrix = *(dKinFitParticles[loc_i]->dCovarianceMatrix); //need to update the values, so don't call Get_CovarianceMatrix() function (returns const matrix) //saves memory this way

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
		if((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex < 0))
		{
			//enclosed decaying particle: the momentum is derived from the momentum of its decay products
				//calculate the uncertainties on this momentum here

		}

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

