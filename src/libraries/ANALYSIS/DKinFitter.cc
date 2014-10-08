#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DKinFitter.h"

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
		//THIS IS CURRENTLY DISABLED
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
	dKinFitStatus = d_KinFitSuccessful;
	dConvergenceChiSqDiff = 0.001;
	dConvergenceChiSqDiff_LastResort = 0.005;

	dMaxKinFitParticlePoolSize = 100;
	dMaxKinFitConstraintVertexPoolSize = 25;
	dMaxKinFitConstraintSpacetimePoolSize = 25;
	dMaxKinFitConstraintP4PoolSize = 25;
	dMaxMatrixDSymPoolSize = 100;
	dMaxLargeMatrixDSymPoolSize = 100;

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

	for(size_t loc_i = 0; loc_i < dLargeMatrixDSymPool_All.size(); ++loc_i)
		delete dLargeMatrixDSymPool_All[loc_i];
}

void DKinFitter::Preallocate_MatrixMemory(void)
{
	//pre-allocate matrix memory
	for(size_t loc_i = 0; loc_i < dMaxMatrixDSymPoolSize; ++loc_i)
	{
		TMatrixDSym* locMatrix = Get_MatrixDSymResource();	
		locMatrix->ResizeTo(7, 7);
	}

	//pre-allocate large matrix memory
	for(size_t loc_i = 0; loc_i < dMaxLargeMatrixDSymPoolSize; ++loc_i)
	{
		TMatrixDSym* locMatrix = Get_LargeMatrixDSymResource();	
		locMatrix->ResizeTo(100, 100);
	}

	Reset_NewEvent();
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

	if(dLargeMatrixDSymPool_All.size() > dMaxLargeMatrixDSymPoolSize)
	{
		for(size_t loc_i = dMaxLargeMatrixDSymPoolSize; loc_i < dLargeMatrixDSymPool_All.size(); ++loc_i)
			delete dLargeMatrixDSymPool_All[loc_i];
		dLargeMatrixDSymPool_All.resize(dMaxLargeMatrixDSymPoolSize);
	}
	dLargeMatrixDSymPool_Available = dLargeMatrixDSymPool_All;
}

void DKinFitter::Reset_NewFit(void)
{
	dKinFitStatus = d_KinFitSuccessful;

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

	dV = Get_LargeMatrixDSymResource();
	dVXi = Get_LargeMatrixDSymResource();
	dVEta = Get_LargeMatrixDSymResource();
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

TMatrixDSym* DKinFitter::Get_LargeMatrixDSymResource(void)
{
	TMatrixDSym* locMatrixDSym;
	if(dLargeMatrixDSymPool_Available.empty())
	{
		locMatrixDSym = new TMatrixDSym();
		dLargeMatrixDSymPool_All.push_back(locMatrixDSym);
	}
	else
	{
		locMatrixDSym = dLargeMatrixDSymPool_Available.back();
		dLargeMatrixDSymPool_Available.pop_back();
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

DKinFitConstraint_P4* DKinFitter::Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint, bool locCloneParticlesFlag)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_P4* locClonedConstraint = Get_KinFitConstraintP4Resource();
	*locClonedConstraint = *locConstraint;

	if(!locCloneParticlesFlag)
		return locClonedConstraint;

	//get or clone member particles
	for(size_t loc_i = 0; loc_i < locConstraint->dInitialParticles.size(); ++loc_i)
		locClonedConstraint->dInitialParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dInitialParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dFinalParticles.size(); ++loc_i)
		locClonedConstraint->dFinalParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dFinalParticles[loc_i]);

	return locClonedConstraint;
}

DKinFitConstraint_Vertex* DKinFitter::Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint, bool locCloneParticlesFlag)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_Vertex* locClonedConstraint = Get_KinFitConstraintVertexResource();
	*locClonedConstraint = *locConstraint;

	if(!locCloneParticlesFlag)
		return locClonedConstraint;

	//get or clone member particles
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
		locClonedConstraint->dFullConstrainParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dFullConstrainParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
		locClonedConstraint->dNoConstrainParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dNoConstrainParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
		locClonedConstraint->dDecayingParticles[loc_i].first = GetOrCreate_ClonedParticle(locConstraint->dDecayingParticles[loc_i].first);

	locClonedConstraint->dDecayingParticlesToAssign.clear();
	set<DKinFitParticle*>::iterator locIterator = locConstraint->dDecayingParticlesToAssign.begin();
	for(; locIterator != locConstraint->dDecayingParticlesToAssign.end(); ++locIterator)
		locClonedConstraint->dDecayingParticlesToAssign.insert(GetOrCreate_ClonedParticle(*locIterator));

	return locClonedConstraint;
}

DKinFitConstraint_Spacetime* DKinFitter::Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint, bool locCloneParticlesFlag)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_Spacetime* locClonedConstraint = Get_KinFitConstraintSpacetimeResource();
	*locClonedConstraint = *locConstraint;

	if(!locCloneParticlesFlag)
		return locClonedConstraint;

	//get or clone member particles
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
		locClonedConstraint->dFullConstrainParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dFullConstrainParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dOnlyConstrainTimeParticles.size(); ++loc_i)
		locClonedConstraint->dOnlyConstrainTimeParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dOnlyConstrainTimeParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
		locClonedConstraint->dNoConstrainParticles[loc_i] = GetOrCreate_ClonedParticle(locConstraint->dNoConstrainParticles[loc_i]);
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
		locClonedConstraint->dDecayingParticles[loc_i].first = GetOrCreate_ClonedParticle(locConstraint->dDecayingParticles[loc_i].first);

	locClonedConstraint->dDecayingParticlesToAssign.clear();
	set<DKinFitParticle*>::iterator locIterator = locConstraint->dDecayingParticlesToAssign.begin();
	for(; locIterator != locConstraint->dDecayingParticlesToAssign.end(); ++locIterator)
		locClonedConstraint->dDecayingParticlesToAssign.insert(GetOrCreate_ClonedParticle(*locIterator));

	return locClonedConstraint;
}

const DKinFitParticle* DKinFitter::Make_DecayingParticle(int locPID, int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_DecayingParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Decaying particle set. Pointer, ID, Q, Mass, Pointer = " << locKinFitParticle << ", " << locPID << ", " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_MissingParticle(int locPID, int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_MissingParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Missing particle set. Pointer, ID, Q, Mass = " << locKinFitParticle << ", " << locPID << ", " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_BeamParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Beam particle set. Pointer, ID, Q, Mass, P3, V3, T = " << locKinFitParticle << ", " << locPID << ", " << locCharge << ", " << locMass << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_TargetParticle(int locPID, int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_TargetParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Target particle set. Pointer, ID, Q, Mass = " << locKinFitParticle << ", " << locPID << ", " << locCharge << ", " << locMass << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Detected particle set. Pointer, ID, Q, Mass, P3, V3, T, pz uncert = " << locKinFitParticle << ", " << locPID << ", " << locCharge << ", " << locMass << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << ", " << sqrt((*locCovarianceMatrix)(3, 3)) << endl;

	return locKinFitParticle;
}

const DKinFitParticle* DKinFitter::Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 5) || (locCovarianceMatrix->GetNcols() != 5))
		return NULL; //is not 5x5

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(0);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_IsNeutralShowerFlag(true);
	locKinFitParticle->Set_Position(locSpacetimeVertex.Vect());
	locKinFitParticle->Set_Time(locSpacetimeVertex.T());
	locKinFitParticle->Set_ShowerEnergy(locShowerEnergy);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
		cout << "DKinFitter: Detected shower set. Pointer, ID, Q, Mass, E, V3, T = 0, " << locKinFitParticle << ", " << locPID << ", " << locMass << ", " << locShowerEnergy << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;

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

DKinFitConstraint_Vertex* DKinFitter::Make_VertexConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, TVector3 locVertexGuess)
{
	deque<DKinFitParticle*> locFullConstrainParticles; //charged particles, decaying particles, beam particles
	deque<pair<DKinFitParticle*, bool> > locDecayingParticles; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
	deque<DKinFitParticle*> locNoConstrainParticles; //missing particles & neutral showers //not used to constrain vertex or time, but fit vertex is set for this particle

	//require all of the tracks to pass through a common point
	//decaying particles are only used to constrain the fit if they are included in exactly two vertex constraints (once as an initial particle, once as a final particle)
		//else they are treated as dNoConstrainParticles
	//will set decaying particles in dFullConstrainParticles and dNoConstrainParticles when ready, but not yet!

	//first check to make sure the inputs are ok //can't tell if enough particles until Resolve_DecayingParticleSpacetimeLinks() //due to decaying particles in 2 constraints
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		if(locInitialParticles[loc_i] == NULL)
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

	//sort particles by how they'll be used by the constraint
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = const_cast<DKinFitParticle*>(locInitialParticles[loc_i]);
		if(locParticle->Get_KinFitParticleType() == d_TargetParticle)
			locNoConstrainParticles.push_back(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locParticle, false));
			else
				locNoConstrainParticles.push_back(locParticle);
		}
		else if(locParticle->Get_KinFitParticleType() == d_BeamParticle)
		{
			//only add if both vx & vy uncertainties are non-zero (else constraints are bad!!)
			const TMatrixDSym& locCovarianceMatrix = *(locParticle->Get_CovarianceMatrix());
			if((locCovarianceMatrix(3, 3) > 0.0) && (locCovarianceMatrix(4, 4) > 0.0)) //include beamline in vertex fit!
				locFullConstrainParticles.push_back(locParticle);
			else
				locNoConstrainParticles.push_back(locParticle);
		}
	}
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = const_cast<DKinFitParticle*>(locFinalParticles[loc_i]);
		if(locParticle->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locParticle, true));
			else
				locNoConstrainParticles.push_back(locParticle);
		}
		else if(locParticle->Get_Charge() == 0)
			locNoConstrainParticles.push_back(locParticle);
		else
			locFullConstrainParticles.push_back(locParticle);
	}

	//create the constraint and set its members
	DKinFitConstraint_Vertex* locKinFitConstraint = Get_KinFitConstraintVertexResource();
	locKinFitConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locKinFitConstraint->Set_DecayingParticles(locDecayingParticles);
	locKinFitConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locKinFitConstraint->Set_CommonVertex(locVertexGuess);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Vertex constraint created. Constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locFullConstrainParticles.size(); ++loc_i)
	 		cout << locFullConstrainParticles[loc_i]->Get_PID() << ", " << locFullConstrainParticles[loc_i]->Get_Charge() << ", " << locFullConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	 		cout << locNoConstrainParticles[loc_i]->Get_PID() << ", " << locNoConstrainParticles[loc_i]->Get_Charge() << ", " << locNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	 		cout << locDecayingParticles[loc_i].first->Get_PID() << ", " << locDecayingParticles[loc_i].first->Get_Charge() << ", " << locDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	return locKinFitConstraint;
}

DKinFitConstraint_Spacetime* DKinFitter::Make_SpacetimeConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locUseRFTimeFlag, TVector3 locVertexGuess, double locCommonTimeGuess)
{
	cout << "ERROR: SPACETIME CONSTRAINTS ARE NOT SUPPORTED YET. RETURNING." << endl;
	return NULL;

	deque<DKinFitParticle*> locFullConstrainParticles; //charged particles, decaying particles, beam particles
	deque<DKinFitParticle*> locOnlyConstrainTimeParticles; //neutral showers //not used to constrain vertex, but fit vertex is used for time constraint
	deque<pair<DKinFitParticle*, bool> > locDecayingParticles; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
	deque<DKinFitParticle*> locNoConstrainParticles; //missing particles //not used to constrain vertex or time, but fit vertex & time are set for this particle

	//require all of the tracks to pass through a common point at a common time
	//decaying particles are only used to constrain the fit if they are included in exactly two vertex constraints (once as an initial particle, once as a final particle)
		//else they are treated as dNoConstrainParticles
	//will set decaying particles in dFullConstrainParticles and dNoConstrainParticles when ready, but not yet!

	//first check to make sure the inputs are ok //can't tell if enough particles until Resolve_DecayingParticleSpacetimeLinks() //due to decaying particles in 2 constraints
	DKinFitParticle* locBeamParticle = NULL;
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		if(locInitialParticles[loc_i] == NULL)
			return NULL;
		else if(locInitialParticles[loc_i]->Get_KinFitParticleType() == d_BeamParticle)
			locBeamParticle = const_cast<DKinFitParticle*>(locInitialParticles[loc_i]);
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

	//sort particles by how they'll be used by the constraint
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = const_cast<DKinFitParticle*>(locInitialParticles[loc_i]);
		if(locParticle->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_TargetParticle)
			locNoConstrainParticles.push_back(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locParticle, false));
			else
				locNoConstrainParticles.push_back(locParticle);
		}
		else if(locParticle->Get_KinFitParticleType() == d_BeamParticle)
		{
			//only add if both vx & vy uncertainties are non-zero (else constraints are bad!!)
			const TMatrixDSym& locCovarianceMatrix = *(locParticle->Get_CovarianceMatrix());
			if((locCovarianceMatrix(3, 3) > 0.0) && (locCovarianceMatrix(4, 4) > 0.0))
				locFullConstrainParticles.push_back(locParticle); //include beamline in vertex fit!
			else
				locNoConstrainParticles.push_back(locParticle);
		}
	}
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = const_cast<DKinFitParticle*>(locFinalParticles[loc_i]);
		if(locParticle->Get_KinFitParticleType() == d_MissingParticle)
			locNoConstrainParticles.push_back(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.push_back(pair<DKinFitParticle*, bool>(locParticle, true));
			else
				locNoConstrainParticles.push_back(locParticle);
		}
		else if(locParticle->Get_Charge() == 0)
			locOnlyConstrainTimeParticles.push_back(locParticle);
		else
			locFullConstrainParticles.push_back(locParticle);
	}

	//create the constraint and set its members
	DKinFitConstraint_Spacetime* locKinFitConstraint = Get_KinFitConstraintSpacetimeResource();
	locKinFitConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locKinFitConstraint->Set_OnlyConstrainTimeParticles(locOnlyConstrainTimeParticles);
	locKinFitConstraint->Set_DecayingParticles(locDecayingParticles);
	locKinFitConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locKinFitConstraint->Set_CommonVertex(locVertexGuess);
	locKinFitConstraint->Set_CommonTime(locCommonTimeGuess);
	locKinFitConstraint->Set_UseRFTimeFlag(locUseRFTimeFlag);
	locKinFitConstraint->dBeamParticle = locBeamParticle;

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Spacetime constraint created. Vertex/Time constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locFullConstrainParticles.size(); ++loc_i)
	 		cout << locFullConstrainParticles[loc_i]->Get_PID() << ", " << locFullConstrainParticles[loc_i]->Get_Charge() << ", " << locFullConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Time-only constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locOnlyConstrainTimeParticles.size(); ++loc_i)
	 		cout << locOnlyConstrainTimeParticles[loc_i]->Get_PID() << ", " << locOnlyConstrainTimeParticles[loc_i]->Get_Charge() << ", " << locOnlyConstrainTimeParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locNoConstrainParticles.size(); ++loc_i)
	 		cout << locNoConstrainParticles[loc_i]->Get_PID() << ", " << locNoConstrainParticles[loc_i]->Get_Charge() << ", " << locNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	 		cout << locDecayingParticles[loc_i].first->Get_PID() << ", " << locDecayingParticles[loc_i].first->Get_Charge() << ", " << locDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	return locKinFitConstraint;
}

DKinFitConstraint_P4* DKinFitter::Make_P4Constraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locConstrainInitialParticleMassFlag)
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

	deque<DKinFitParticle*> locNonConstInitialParticles, locNonConstFinalParticles;
	for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
		locNonConstInitialParticles.push_back(const_cast<DKinFitParticle*>(locInitialParticles[loc_i]));
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
		locNonConstFinalParticles.push_back(const_cast<DKinFitParticle*>(locFinalParticles[loc_i]));

	//create the constraint and set its members
	DKinFitConstraint_P4* locKinFitConstraint = Get_KinFitConstraintP4Resource();
	locKinFitConstraint->Set_InitialParticles(locNonConstInitialParticles);
	locKinFitConstraint->Set_FinalParticles(locNonConstFinalParticles);
	locKinFitConstraint->Set_ConstrainInitialParticleMassFlag(locConstrainInitialParticleMassFlag);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: P4 constraint created. Initial-state particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locInitialParticles.size(); ++loc_i)
	 		cout << locInitialParticles[loc_i]->Get_PID() << ", " << locInitialParticles[loc_i]->Get_Charge() << ", " << locInitialParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Final-state particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	 		cout << locFinalParticles[loc_i]->Get_PID() << ", " << locFinalParticles[loc_i]->Get_Charge() << ", " << locFinalParticles[loc_i]->Get_Mass() << endl;
	}

	return locKinFitConstraint;
}

const DKinFitConstraint_P4* DKinFitter::Set_Constraint(const DKinFitConstraint_P4* locConstraint)
{
	//register a constraint for upcoming fit

	//clone particles & constraints so that they can be modified during the fit without overwriting the originals
	DKinFitConstraint_P4* locClonedConstraint = Clone_KinFitConstraint_P4(locConstraint, true);
	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locClonedConstraint));

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: P4 constraint set. Initial-state particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dInitialParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dInitialParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dInitialParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dInitialParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Final-state particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dFinalParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dFinalParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dFinalParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dFinalParticles[loc_i]->Get_Mass() << endl;
	}

	return locClonedConstraint;
}

const DKinFitConstraint_Vertex* DKinFitter::Set_Constraint(const DKinFitConstraint_Vertex* locConstraint)
{
	//register a constraint for upcoming fit

	//clone particles & constraints so that they can be modified during the fit without overwriting the originals
	DKinFitConstraint_Vertex* locClonedConstraint = Clone_KinFitConstraint_Vertex(locConstraint, true);
	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locClonedConstraint));

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Vertex constraint set. Constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dFullConstrainParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dNoConstrainParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dDecayingParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dDecayingParticles[loc_i].first->Get_PID() << ", " << locClonedConstraint->dDecayingParticles[loc_i].first->Get_Charge() << ", " << locClonedConstraint->dDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	return locClonedConstraint;
}

const DKinFitConstraint_Spacetime* DKinFitter::Set_Constraint(const DKinFitConstraint_Spacetime* locConstraint)
{
	//register a constraint for upcoming fit
	cout << "ERROR: SPACETIME CONSTRAINTS ARE NOT SUPPORTED YET. RETURNING." << endl;
	return NULL;

	//clone particles & constraints so that they can be modified during the fit without overwriting the originals
	DKinFitConstraint_Spacetime* locClonedConstraint = Clone_KinFitConstraint_Spacetime(locConstraint, true);
	dKinFitConstraints.push_back(static_cast<DKinFitConstraint*>(locClonedConstraint));

	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Spacetime constraint set. Vertex/Time constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dFullConstrainParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dFullConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Time-only constrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dOnlyConstrainTimeParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dOnlyConstrainTimeParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dOnlyConstrainTimeParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dOnlyConstrainTimeParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Unconstrained particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dNoConstrainParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_PID() << ", " << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_Charge() << ", " << locClonedConstraint->dNoConstrainParticles[loc_i]->Get_Mass() << endl;
		cout << "DKinFitter: Decaying particle PID's, q's, masses: " << endl;
		for(size_t loc_i = 0; loc_i < locClonedConstraint->dDecayingParticles.size(); ++loc_i)
	 		cout << locClonedConstraint->dDecayingParticles[loc_i].first->Get_PID() << ", " << locClonedConstraint->dDecayingParticles[loc_i].first->Get_Charge() << ", " << locClonedConstraint->dDecayingParticles[loc_i].first->Get_Mass() << endl;
	}

	return locClonedConstraint;
}

bool DKinFitter::Prepare_Constraint(DKinFitConstraint_P4* locConstraint) const
{
	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locConstraint->dInitialParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dInitialParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		size_t locNumP4Constraints = locParticle->Get_NumP4Constraints();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one P4 constraint.  Constraint not added." << endl;
			return false;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two P4 constraints.  Constraint not added." << endl;
			return false;
		}
	}
	for(size_t loc_i = 0; loc_i < locConstraint->dFinalParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dFinalParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		size_t locNumP4Constraints = locParticle->Get_NumP4Constraints();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumP4Constraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one P4 constraint.  Constraint not added." << endl;
			return false;
		}
		else if(locNumP4Constraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two P4 constraints.  Constraint not added." << endl;
			return false;
		}
	}

	//mark constraint in particles
	for(size_t loc_i = 0; loc_i < locConstraint->dInitialParticles.size(); ++loc_i)
		locConstraint->dInitialParticles[loc_i]->Add_P4Constraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dFinalParticles.size(); ++loc_i)
		locConstraint->dFinalParticles[loc_i]->Add_P4Constraint(locConstraint);

	return true;
}

bool DKinFitter::Prepare_Constraint(DKinFitConstraint_Vertex* locConstraint) const
{
	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dFullConstrainParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		size_t locNumConstraints = locParticle->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumConstraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one vertex constraint.  Constraint not added." << endl;
			return false;
		}
		else if(locNumConstraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two vertex constraints.  Constraint not added." << endl;
			return false;
		}
	}
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dDecayingParticles[loc_i].first;
		size_t locNumConstraints = locParticle->Get_NumVertexFits();
		if(locNumConstraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two vertex constraints.  Constraint not added." << endl;
			return false;
		}
	}

	//set vertex constraint flags
	TVector3 locMomentum;
	unsigned short int locVertexConstraintFlag;
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = locConstraint->dFullConstrainParticles[loc_i];
		locMomentum = locParticle->Get_Momentum();
		if(fabs(locMomentum.Pz()) > fabs(locMomentum.Px()))
			locVertexConstraintFlag = (fabs(locMomentum.Pz()) > fabs(locMomentum.Py())) ? 1 : 2;
		else
			locVertexConstraintFlag = (fabs(locMomentum.Px()) > fabs(locMomentum.Py())) ? 3 : 2;
		locParticle->Set_VertexConstraintFlag(locVertexConstraintFlag);
	}

	//set momentum of neutral showers that have 5x5 covariance matrix (shower energy input instead of p3)
	for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = locConstraint->dNoConstrainParticles[loc_i];
		if(!locParticle->Get_IsNeutralShowerFlag())
			continue; //only do for neutral showers

		double locE = locParticle->Get_ShowerEnergy();
		double locMass = locParticle->Get_Mass();
		double locPMag = sqrt(locE*locE - locMass*locMass);
		TVector3 locMomentum = locParticle->Get_Position() - locParticle->Get_CommonVertex();
		locMomentum.SetMag(locPMag);
		locParticle->Set_Momentum(locMomentum);
	}

	//add constraint to particles
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
		locConstraint->dFullConstrainParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
		locConstraint->dNoConstrainParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
		locConstraint->dDecayingParticles[loc_i].first->Add_CommonVertexAndOrTimeConstraint(locConstraint);

	return true;
}

bool DKinFitter::Prepare_Constraint(DKinFitConstraint_Spacetime* locConstraint) const
{
	//prepare a constraint for upcoming fit
	cout << "ERROR: SPACETIME CONSTRAINTS ARE NOT SUPPORTED YET. RETURNING." << endl;
	return false;

	//enforce maximum #constraints per particle: 1 for non-decaying, 2 for decaying
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dFullConstrainParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		size_t locNumConstraints = locParticle->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumConstraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one spacetime constraint.  Constraint not added." << endl;
			return false;
		}
		else if(locNumConstraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two spacetime constraints.  Constraint not added." << endl;
			return false;
		}
	}
	for(size_t loc_i = 0; loc_i < locConstraint->dOnlyConstrainTimeParticles.size(); ++loc_i) //neutral showers
	{
		const DKinFitParticle* locParticle = locConstraint->dOnlyConstrainTimeParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();
		size_t locNumConstraints = locParticle->Get_NumVertexFits();
		if((locKinFitParticleType != d_DecayingParticle) && (locNumConstraints > 0))
		{
			cout << "ERROR: Non-decaying particle cannot be used in more than one spacetime constraint.  Constraint not added." << endl;
			return false;
		}
	}
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
	{
		const DKinFitParticle* locParticle = locConstraint->dDecayingParticles[loc_i].first;
		size_t locNumConstraints = locParticle->Get_NumVertexFits();
		if(locNumConstraints > 1)
		{
			cout << "ERROR: Decaying particle cannot be used in more than two spacetime constraints.  Constraint not added." << endl;
			return false;
		}
	}

	//set momentum of neutral showers that have 5x5 covariance matrix (shower energy input instead of p3)
	for(size_t loc_i = 0; loc_i < locConstraint->dOnlyConstrainTimeParticles.size(); ++loc_i)
	{
		DKinFitParticle* locParticle = locConstraint->dOnlyConstrainTimeParticles[loc_i];
		double locE = locParticle->Get_ShowerEnergy();
		double locMass = locParticle->Get_Mass();
		double locPMag = sqrt(locE*locE - locMass*locMass);
		TVector3 locMomentum = locParticle->Get_Position() - locParticle->Get_CommonVertex();
		locMomentum.SetMag(locPMag);
		locParticle->Set_Momentum(locMomentum);
	}

	//add constraint to particles
	for(size_t loc_i = 0; loc_i < locConstraint->dFullConstrainParticles.size(); ++loc_i)
		locConstraint->dFullConstrainParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dOnlyConstrainTimeParticles.size(); ++loc_i)
		locConstraint->dOnlyConstrainTimeParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
		locConstraint->dNoConstrainParticles[loc_i]->Add_CommonVertexAndOrTimeConstraint(locConstraint);
	for(size_t loc_i = 0; loc_i < locConstraint->dDecayingParticles.size(); ++loc_i)
		locConstraint->dDecayingParticles[loc_i].first->Add_CommonVertexAndOrTimeConstraint(locConstraint);

	return true;
}

bool DKinFitter::Sort_Constraints(const deque<DKinFitConstraint*>& locOriginalConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints)
{
#ifdef VTRACE
	VT_TRACER("DKinFitter::Sort_Constraints()");
#endif
	if(dDebugLevel > 10)
		cout << "DKinFitter: Sort constraints: Clone constraints." << endl;

	//clone constraints & particles
	deque<DKinFitConstraint*> locClonedConstraints;
	map<DKinFitConstraint_VertexBase*, DKinFitConstraint_VertexBase*> locCloneToOriginalVertexConstraintMap;
	map<DKinFitConstraint_P4*, DKinFitConstraint_P4*> locCloneToOriginalP4ConstraintMap;
	for(size_t loc_i = 0; loc_i < locOriginalConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locOriginalConstraints[loc_i]);
		if(locP4Constraint != NULL)
		{
			locClonedConstraints.push_back(dynamic_cast<DKinFitConstraint*>(Clone_KinFitConstraint_P4(locP4Constraint, true)));
			locCloneToOriginalP4ConstraintMap[dynamic_cast<DKinFitConstraint_P4*>(locClonedConstraints.back())] = locP4Constraint;
			continue;
		}
		DKinFitConstraint_Vertex* locVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(locOriginalConstraints[loc_i]);
		if(locVertexConstraint != NULL)
		{
			locClonedConstraints.push_back(dynamic_cast<DKinFitConstraint*>(Clone_KinFitConstraint_Vertex(locVertexConstraint, true)));
			locCloneToOriginalVertexConstraintMap[dynamic_cast<DKinFitConstraint_VertexBase*>(locClonedConstraints.back())] = dynamic_cast<DKinFitConstraint_VertexBase*>(locVertexConstraint);
			continue;
		}
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locOriginalConstraints[loc_i]);
		if(locSpacetimeConstraint != NULL)
		{
			locClonedConstraints.push_back(dynamic_cast<DKinFitConstraint*>(Clone_KinFitConstraint_Spacetime(locSpacetimeConstraint, true)));
			locCloneToOriginalVertexConstraintMap[dynamic_cast<DKinFitConstraint_VertexBase*>(locClonedConstraints.back())] = dynamic_cast<DKinFitConstraint_VertexBase*>(locSpacetimeConstraint);
			continue;
		}
	}

	if(dDebugLevel > 10)
		cout << "DKinFitter: Sort constraints: Resolve constraints." << endl;

	deque<DKinFitConstraint_VertexBase*> locSortedVertexConstraints;
	if(!Resolve_Constraints(locClonedConstraints, locSortedVertexConstraints, true))
		return false;

	if(dDebugLevel > 10)
		cout << "DKinFitter: Sort constraints: Group constraints." << endl;

	if(!Group_Constraints(locSortedVertexConstraints, locSortedConstraints))
		return false;

	//loop over locSortedConstraints, replacing cloned contents with originals
	for(size_t loc_i = 0; loc_i < locSortedConstraints.size(); ++loc_i)
	{
		//vertex
		DKinFitConstraint_VertexBase* locClonedVertexConstraint = locSortedConstraints[loc_i].first;
		locSortedConstraints[loc_i].first = locCloneToOriginalVertexConstraintMap[locClonedVertexConstraint];

		//p4
		set<DKinFitConstraint_P4*> locClonedP4Constraints = locSortedConstraints[loc_i].second;
		locSortedConstraints[loc_i].second.clear();
		set<DKinFitConstraint_P4*>::iterator locIterator = locClonedP4Constraints.begin();
		for(; locIterator != locClonedP4Constraints.end(); ++locIterator)
		{
			DKinFitConstraint_P4* locClonedP4Constraint = *locIterator;
			locSortedConstraints[loc_i].second.insert(locCloneToOriginalP4ConstraintMap[locClonedP4Constraint]);
		}
	}

	//loop over locCloneToOriginalP4ConstraintMap, modifying originals with updated constrain-mass flags
	map<DKinFitConstraint_P4*, DKinFitConstraint_P4*>::iterator locMapIterator = locCloneToOriginalP4ConstraintMap.begin();
	for(; locMapIterator != locCloneToOriginalP4ConstraintMap.end(); ++locMapIterator)
		locMapIterator->second->Set_ConstrainInitialParticleMassFlag(locMapIterator->first->Get_ConstrainInitialParticleMassFlag());

	return true;
}

bool DKinFitter::Resolve_Constraints(void)
{
	deque<DKinFitConstraint_VertexBase*> locSortedVertexConstraints;
	return Resolve_Constraints(dKinFitConstraints, locSortedVertexConstraints, false);
}

bool DKinFitter::Resolve_Constraints(const deque<DKinFitConstraint*>& locConstraints, deque<DKinFitConstraint_VertexBase*>& locSortedVertexConstraints, bool locSortOnlyFlag) const
{
#ifdef VTRACE
	VT_TRACER("DKinFitter::Resolve_Constraints()");
#endif
	if(dDebugLevel > 10)
		cout << "DKinFitter: Resolve constraints: Prepare " << locConstraints.size() << " constraints." << endl;

	//prepare constraints
	for(size_t loc_i = 0; loc_i < locConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locConstraints[loc_i]);
		if(locP4Constraint != NULL)
		{
			if(!Prepare_Constraint(locP4Constraint))
				return false;
			continue;
		}

		DKinFitConstraint_Vertex* locVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(locConstraints[loc_i]);
		if(locVertexConstraint != NULL)
		{
			if(!Prepare_Constraint(locVertexConstraint))
				return false;
			continue;
		}
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locConstraints[loc_i]);
		if(locSpacetimeConstraint != NULL)
		{
			if(!Prepare_Constraint(locSpacetimeConstraint))
				return false;
			continue;
		}
	}

	if(dDebugLevel > 10)
		cout << "DKinFitter: Resolve constraints: Constraints Prepared." << endl;

	locSortedVertexConstraints.clear();
	if(!Resolve_DecayingParticleSpacetimeLinks(locConstraints, locSortedVertexConstraints, locSortOnlyFlag))
		return false;

	if(dDebugLevel > 10)
		cout << "DKinFitter: Resolve constraints: Spacetime resolved." << endl;

	if(!Resolve_P4Constraints(locConstraints, locSortOnlyFlag))
		return false;

	if(dDebugLevel > 10)
		cout << "DKinFitter: Resolve constraints: P4 resolved." << endl;

	if(!Resolve_P4MassConstraints(locConstraints, locSortOnlyFlag))
		return false;

	if(dDebugLevel > 10)
		cout << "DKinFitter: Resolve constraints: P4/Mass resolved." << endl;

	return true;
}

bool DKinFitter::Resolve_DecayingParticleSpacetimeLinks(const deque<DKinFitConstraint*>& locKinFitConstraints, deque<DKinFitConstraint_VertexBase*>& locSortedConstraints, bool locSortOnlyFlag) const
{
	//resolve links between vertex & time fits (decaying particles), and return sorted constraints

	//if locSortOnlyFlag = true, then cut-out invalid vertex constraints, rather than returning false (done if false)
	locSortedConstraints.clear();

	//build deque of DKinFitConstraint_VertexBase to sort through
	deque<DKinFitConstraint_VertexBase*> locVertexConstraintsToSort;
	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_VertexBase* locKinFitConstraint_VertexBase = dynamic_cast<DKinFitConstraint_VertexBase*>(locKinFitConstraints[loc_i]);
		if(locKinFitConstraint_VertexBase != NULL)
			locVertexConstraintsToSort.push_back(locKinFitConstraint_VertexBase);
	}

	//loop over vertex-constraints-to-sort:
		//find which constraints decaying particles should be defined-by/constrained-to
		//find order in which constraints need to be constrained
	deque<DKinFitConstraint_VertexBase*>::iterator locSortIterator = locVertexConstraintsToSort.begin();
	bool locProgessMadeFlag = false;
	bool locLastResortFlag = false;
	bool locInitPassFlag = true;
	while(true)
	{
		if(locVertexConstraintsToSort.empty())
			break; //all vertex constraints setup successfully
		if(locSortIterator == locVertexConstraintsToSort.end())
		{
			if((!locProgessMadeFlag) && (!locInitPassFlag))
			{
				if(locLastResortFlag)
				{
					if(!locSortOnlyFlag)
						cout << "ERROR: NOT ENOUGH PARTICLES TO CONSTRAIN VERTEX." << endl;
					//no progress made, and was on last resort: cannot constrain remaining vertices
					return locSortOnlyFlag; //this is ok if only sorting (skim remaining constraints), but bad if fit time
				}
				locLastResortFlag = true; //no progress, now on last resort
			}
			locInitPassFlag = false;
			locSortIterator = locVertexConstraintsToSort.begin();
			locProgessMadeFlag = false;
			continue;
		}

		DKinFitConstraint_VertexBase* locConstraint = *locSortIterator;
		if(locConstraint->dFullConstrainParticles.size() < 2)
		{
			++locSortIterator;
			continue;
		}

		//for the init pass through all of the constraints, only accept constraints that don't have 2 decaying/missing particles
			//this way, any decaying particles at these vertices have p4's that are defined locally (instead of across many steps)
		if(locInitPassFlag)
		{
			size_t locNumMissingDecayingParticles = locConstraint->Get_DecayingParticles().size();
			for(size_t loc_i = 0; loc_i < locConstraint->dNoConstrainParticles.size(); ++loc_i)
			{
				if(locConstraint->dNoConstrainParticles[loc_i]->Get_KinFitParticleType() != d_MissingParticle)
					continue;
				++locNumMissingDecayingParticles;
				break;
			}
			if(locNumMissingDecayingParticles >= 2)
			{
				++locSortIterator;
				continue;
			}
		}

		if((locConstraint->dDecayingParticlesToAssign.size() >= 2) && (!locLastResortFlag))
		{
			++locSortIterator;
			continue; //constrain decaying particles with other constraints first, if at all possible (unless last resort!)
		}

		//any remaining decaying particles can now be defined here (if not already), and added as constraints to their other vertex fits
		set<DKinFitParticle*>::iterator locIterator = locConstraint->dDecayingParticlesToAssign.begin();
		for(; locIterator != locConstraint->dDecayingParticlesToAssign.end(); ++locIterator)
		{
			DKinFitParticle* locParticle = *locIterator;
			if(locParticle->Get_DefinedAtVertexAndOrTimeConstraint() != NULL)
				continue; //already defined, meaning already added as a constraint previously as well

			//set particle info
			locParticle->Set_Position(locConstraint->Get_CommonVertex());
			DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locConstraint);
			if(locSpacetimeConstraint != NULL)
				locParticle->Set_Time(locSpacetimeConstraint->Get_CommonTime());

			//set whether it's defined at it the production vertex or decay vertex
				//note: if a decaying particle is not in a vertex fit, then this quantity doesn't matter
			bool locProductionVertexFlag = true; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
			deque<pair<DKinFitParticle*, bool> > locDecayingParticles = locConstraint->dDecayingParticles;
			for(size_t loc_k = 0; loc_k < locDecayingParticles.size(); ++loc_k)
			{
				if(locDecayingParticles[loc_k].first != locParticle)
					continue;
				locProductionVertexFlag = locDecayingParticles[loc_k].second;
				break;
			}
			locParticle->Set_DecayingParticleAtProductionVertexFlag(locProductionVertexFlag);

			//define it, add to constraints of other vertex fits
			locConstraint->Add_NoConstrainParticle(locParticle);
			deque<DKinFitConstraint_VertexBase*> locNewVertexConstraints = locParticle->dCommonVertexAndOrTimeConstraints;
			for(size_t loc_k = 0; loc_k < locNewVertexConstraints.size(); ++loc_k)
			{
				if(locNewVertexConstraints[loc_k] == locConstraint)
					continue;
				//in multiple vertex fits
				locNewVertexConstraints[loc_k]->Add_FullConstrainParticle(locParticle);
				locNewVertexConstraints[loc_k]->dDecayingParticlesToAssign.erase(locParticle);
				//make sure also in a p4 fit so its p4 is defined, else cannot constrain (if linked)
				if((!locParticle->Get_IsInP4FitFlag()) && dLinkVerticesFlag)
				{
					//unrecoverable, even if only-sorting: problem is with the p4 constraints
					cout << "ERROR in DKinFitter: Decaying particle constrained in a vertex or spacetime fit with unconstrained momentum!!  Exiting." << endl;
					return false; //decaying but constrained in a vertex or spacetime fit with unconstrained momentum
				}
			}
		}
		locConstraint->dDecayingParticlesToAssign.clear();
		locSortedConstraints.push_back(*locSortIterator);
		locProgessMadeFlag = true;
		locLastResortFlag = false;
		locSortIterator = locVertexConstraintsToSort.erase(locSortIterator);
	}

	return true;
}

bool DKinFitter::Resolve_P4Constraints(const deque<DKinFitConstraint*>& locKinFitConstraints, bool locSortOnlyFlag) const
{
	deque<DKinFitConstraint_P4*> locP4Constraints;
	set<DKinFitParticle*> locParticlesInP4Constraints;
	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(locKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 == NULL)
			continue;
		locP4Constraints.push_back(locKinFitConstraint_P4);

		for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dInitialParticles.size(); ++loc_j)
			locParticlesInP4Constraints.insert(locKinFitConstraint_P4->dInitialParticles[loc_j]);
		for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
			locParticlesInP4Constraints.insert(locKinFitConstraint_P4->dFinalParticles[loc_j]);
	}

	//snag unconstrained particles
	deque<DKinFitParticle*> locUnconstrainedParticles;
	set<DKinFitParticle*>::iterator locIterator = locParticlesInP4Constraints.begin();
	for(; locIterator != locParticlesInP4Constraints.end(); ++locIterator)
	{
		DKinFitParticle* locKinFitParticle = *locIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		//make sure there are no neutral showers used in p4 constraints that are not included in a vertex fit
		if(locKinFitParticle->Get_IsNeutralShowerFlag() && (!locKinFitParticle->Get_IsInVertexOrSpacetimeFitFlag()) && (!locSortOnlyFlag))
		{
			cout << "ERROR in DKinFitter: Detected neutral shower in a P4 fit but not in a vertex or spacetime fit: P3 is undefined!!  Exiting." << endl;
			return false; //detected neutral shower in a p4 fit but not in a vertex fit: p3 is undefined!!
		}
		if((locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
			continue;
		locUnconstrainedParticles.push_back(locKinFitParticle);
	}

	//get p4 guesses for each decaying/missing particle
	//loop through the p4 constraints, find constraints which contain only one missing or decaying particle: used as starting point for assigning them to different constraints
	deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> > locConstrainableParticles;
	deque<const DKinFitParticle*> locConstrainedParticles;
	while(locConstrainedParticles.size() < locUnconstrainedParticles.size())
	{
		if(locConstrainableParticles.empty())
		{
			if(!Find_ConstrainableParticles(locP4Constraints, locConstrainableParticles, locConstrainedParticles))
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

			DKinFitParticle* locParticleToConstrain = locConstrainableParticles.back().first;
			DKinFitConstraint_P4* locKinFitConstraint_P4 = locConstrainableParticles.back().second;
			locKinFitConstraint_P4->Set_ConstrainedP4Particle(locParticleToConstrain);
			TLorentzVector locP4;
			Constrain_Particle(locParticleToConstrain, locKinFitConstraint_P4, locP4);
			locParticleToConstrain->Set_Momentum(locP4.Vect());
			if(locParticleToConstrain->Get_PID() == 0)
				locParticleToConstrain->Set_Mass(locP4.M());
			if(dDebugLevel > 5)
			{
				if(locParticleToConstrain->Get_DecayingParticleAtProductionVertexFlag())
					cout << "particle is defined at its production vertex (possibly constrained to its decay vertex)" << endl;
				else
					cout << "particle is defined at its decay vertex (possibly constrained to its production vertex)" << endl;
			}

			//set vertex constraint flag if decaying particle & in 2 vertex fits
			if(locParticleToConstrain->Get_NumVertexFits() == 2)
			{
				unsigned short int locVertexConstraintFlag;
				if(fabs(locP4.Pz()) > fabs(locP4.Px()))
					locVertexConstraintFlag = (fabs(locP4.Pz()) > fabs(locP4.Py())) ? 1 : 2;
				else
					locVertexConstraintFlag = (fabs(locP4.Px()) > fabs(locP4.Py())) ? 3 : 2;
				locParticleToConstrain->Set_VertexConstraintFlag(locVertexConstraintFlag);
			}

			locConstrainedParticles.push_back(locConstrainableParticles.back().first);
			locConstrainableParticles.pop_back();
		}
	}

	return true;
}

bool DKinFitter::Find_ConstrainableParticles(const deque<DKinFitConstraint_P4*>& locP4Constraints, deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles) const
{
	for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
	{
		if(locP4Constraints[loc_i]->dConstrainedP4Particle != NULL)
			continue; //constraint already set
		size_t locNumParticlesNeedToBeConstrained = 0;
		DKinFitParticle* locConstrainableParticle = NULL;
		deque<DKinFitParticle*> locTempParticles = locP4Constraints[loc_i]->dInitialParticles;
		for(size_t loc_j = 0; loc_j < locTempParticles.size(); ++loc_j)
		{
			DKinFitParticleType locKinFitParticleType = locTempParticles[loc_j]->Get_KinFitParticleType();
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
		locTempParticles = locP4Constraints[loc_i]->dFinalParticles;
		for(size_t loc_j = 0; loc_j < locTempParticles.size(); ++loc_j)
		{
			DKinFitParticleType locKinFitParticleType = locTempParticles[loc_j]->Get_KinFitParticleType();
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
			locConstrainableParticles.push_back(pair<DKinFitParticle*, DKinFitConstraint_P4*>(locConstrainableParticle, locP4Constraints[loc_i]));
	}
	return (!locConstrainableParticles.empty());
}

void DKinFitter::Constrain_Particle(DKinFitParticle* locParticleToConstrain, DKinFitConstraint_P4* locConstraint, TLorentzVector& locP4) const
{
	if(dDebugLevel > 5)
		cout << "particle to constrain PID, q, mass = " << locParticleToConstrain->Get_PID() << ", " << locParticleToConstrain->Get_Charge() << ", " << locParticleToConstrain->Get_Mass() << endl;

	locP4.SetXYZT(0.0, 0.0, 0.0, 0.0);
	bool locConstrainedParticleIsInInitialState = false;
	for(size_t loc_j = 0; loc_j < locConstraint->dInitialParticles.size(); ++loc_j)
	{
		DKinFitParticle* locKinFitParticle = locConstraint->dInitialParticles[loc_j];
		if(dDebugLevel > 20)
			cout << "init particle PID, q, mass = " << locKinFitParticle->Get_PID() << ", " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		if(locKinFitParticle == locParticleToConstrain)
			locConstrainedParticleIsInInitialState = true;
		else
			locP4 += locKinFitParticle->Get_P4();
	}
	for(size_t loc_j = 0; loc_j < locConstraint->dFinalParticles.size(); ++loc_j)
	{
		DKinFitParticle* locKinFitParticle = locConstraint->dFinalParticles[loc_j];
		if(dDebugLevel > 20)
			cout << "final particle PID, q, mass = " << locKinFitParticle->Get_PID() << ", " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		if(locKinFitParticle != locParticleToConstrain)
			locP4 -= locKinFitParticle->Get_P4();
	}

	//set p3 guess
	if(locConstrainedParticleIsInInitialState)
		locP4 *= -1.0;

	if(dDebugLevel > 5)
		cout << "particle to constrain: flag, pxyzE = " << locConstrainedParticleIsInInitialState << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	locConstraint->Set_ConstrainedParticleIsInInitialStateFlag(locConstrainedParticleIsInInitialState);
}

bool DKinFitter::Resolve_P4MassConstraints(const deque<DKinFitConstraint*>& locKinFitConstraints, bool locSortOnlyFlag) const
{
	if(!Resolve_InclusiveP4(locKinFitConstraints))
		return false; //weird p4 constraints

	//remove/return-false if mass + p4 constraints too over-constrain the system (will result in noninvertable matrices)

	//first find "the" p4 constraint
	DKinFitConstraint_P4* locTheP4Constraint = NULL;
	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locKinFitConstraints[loc_i]);
		if(locP4Constraint == NULL)
			continue;
		if(!locP4Constraint->Get_IsActualP4ConstraintFlag())
			continue;
		locTheP4Constraint = locP4Constraint;
		break;
	}
	if(locTheP4Constraint == NULL)
		return true; //inclusive fit: not a problem
	locTheP4Constraint->dConstrainMassFlag = false;

	//ok, now loop over the p4 constraint: if it has a particle with non-zero cov matrix entries that is not in a p4 constraint, then we'll be OK
	for(size_t loc_i = 0; loc_i < locTheP4Constraint->dInitialParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locTheP4Constraint->dInitialParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle))
			continue;
		if(locKinFitParticle->Get_CovarianceMatrix() == NULL)
			continue;
		if(!(fabs((*locKinFitParticle->Get_CovarianceMatrix())(2, 2)) > 0.0))
			continue; //error on pz is zero
		return true; //have a particle in the p4 constraint with non-zero errors that is not in the mass constraint
	}

	deque<DKinFitParticle*> locDecayingParticles;
	for(size_t loc_i = 0; loc_i < locTheP4Constraint->dFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locTheP4Constraint->dFinalParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if(locKinFitParticleType == d_DecayingParticle)
		{
			locDecayingParticles.push_back(locKinFitParticle);
			continue;
		}
		if(locKinFitParticleType == d_MissingParticle)
			continue;
		if(locKinFitParticle->Get_CovarianceMatrix() == NULL)
			continue;
		if(!(fabs((*locKinFitParticle->Get_CovarianceMatrix())(2, 2)) > 0.0))
			continue; //error on pz is zero
		return true; //have a particle in the p4 constraint with non-zero errors that is not in the mass constraint
	}

	if(locDecayingParticles.empty())
		return false; //not a single particle with non-zero errors available

	//ok, we might have a problem: make sure that a mass constraint is NOT applied to at least one of the decaying particles
	for(size_t loc_i = 0; loc_i < locDecayingParticles.size(); ++loc_i)
	{
		const DKinFitConstraint_P4* locConstraintAsInitial = locDecayingParticles[loc_i]->Get_P4ConstraintWhenInitial();
		if(locConstraintAsInitial->dConstrainMassFlag)
			continue;
		return true; //a decaying particle does not have a mass constraint applied: we're ok
	}

	//ok, this won't work. if only sorting, disable a mass constraint. else return false
	if(locSortOnlyFlag)
	{
		DKinFitConstraint_P4* locConstraintAsInitial = const_cast<DKinFitConstraint_P4*>(locDecayingParticles[0]->Get_P4ConstraintWhenInitial());
		locConstraintAsInitial->dConstrainMassFlag = false; //could disable any of them, but just choose this one
		if(dDebugLevel > 0)
			cout << "DKinFitter: System too over-constrained, removed mass constraint on PID = " << locDecayingParticles[0]->Get_PID() << endl;
		return true;
	}

	cout << "ERROR: CANNOT APPLY OVERALL P4 CONSTRAINT AND THESE MASS CONSTRAINTS AT THE SAME TIME: SYSTEM IS TOO OVER-CONSTRAINED. RETURNING FALSE" << endl;
	return false;
}

bool DKinFitter::Resolve_InclusiveP4(const deque<DKinFitConstraint*>& locKinFitConstraints) const
{
	//if is an inclusive fit, will mark mass constraints as missing if necessary
	//if not an inclusive fit, will mark which constraint is the actual p4 constraint

	//see if there is a missing particle in a p4 constraint that has PID = 0
	DKinFitConstraint_P4* locInclusiveP4Constraint = NULL;
	bool locUnknownParentFoundFlag = false;
	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locKinFitConstraints[loc_i]);
		if(locP4Constraint == NULL)
			continue;

		for(size_t loc_i = 0; loc_i < locP4Constraint->dInitialParticles.size(); ++loc_i)
		{
			DKinFitParticle* locKinFitParticle = locP4Constraint->dFinalParticles[loc_i];
			if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
				continue;
			if(locKinFitParticle->Get_PID() == 0)
			{
				locUnknownParentFoundFlag = true;
				locInclusiveP4Constraint = locP4Constraint;
				if(dDebugLevel > 20)
					cout << "Decaying particle with unknown PID found: is inclusive-p4 fit." << endl;
				locInclusiveP4Constraint->dConstrainMassByInvariantMassFlag = false;
			}
			break;
		}
		if(locUnknownParentFoundFlag)
			break;

		bool locMissingParticleFoundFlag = false;
		for(size_t loc_i = 0; loc_i < locP4Constraint->dFinalParticles.size(); ++loc_i)
		{
			DKinFitParticle* locKinFitParticle = locP4Constraint->dFinalParticles[loc_i];
			if(locKinFitParticle->Get_KinFitParticleType() != d_MissingParticle)
				continue;
			locMissingParticleFoundFlag = true;
			if(locKinFitParticle->Get_PID() == 0)
			{
				locInclusiveP4Constraint = locP4Constraint;
				//if constraining the mass of this initial particle, must do by missing mass (X is a decay product)
				if(dDebugLevel > 20)
					cout << "Missing particle with unknown PID found: is inclusive-p4 fit." << endl;
				locInclusiveP4Constraint->dConstrainMassByInvariantMassFlag = false;
				if(dDebugLevel > 5)
					cout << "p4 constraint with pid " << locInclusiveP4Constraint->dInitialParticles[0]->Get_PID() << " as parent marked as a missing-mass constraint (if mass constrained at all)" << endl;
			}
			break;
		}
		if(locMissingParticleFoundFlag)
			break;
	}

	bool locAreP4ConstraintsFlag = false;
	for(size_t loc_i = 0; loc_i < locKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(locKinFitConstraints[loc_i]);
		if(locP4Constraint == NULL)
			continue;
		locAreP4ConstraintsFlag = true;
		DKinFitParticle* locKinFitParticle = locP4Constraint->dInitialParticles[0];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_BeamParticle) || (locKinFitParticleType == d_DetectedParticle))
		{
			//initial particle is beam particle or detected (is actually a decaying particle treated as detected) particle: p4 constraint instead of mass constraint (if not inclusive)
			locP4Constraint->dConstrainMassFlag = false;
			if(locInclusiveP4Constraint == NULL)
			{
				if(dDebugLevel > 5)
					cout << "p4 constraint with pid " << locKinFitParticle->Get_PID() << " as parent marked as full-p4 constraint" << endl;
				locP4Constraint->dIsActualP4ConstraintFlag = true;
				return true; //is not inclusive: we're done
			}
			if(dDebugLevel > 5)
				cout << "p4 constraint with pid " << locKinFitParticle->Get_PID() << " as parent marked as no-mass constraint" << endl;
		}
		else if((locKinFitParticleType == d_DecayingParticle) && (locKinFitParticle->Get_NumP4Constraints() == 1))
		{
			//initial particle is open-ended decaying particle: p4 constraint instead of mass constraint (if not inclusive)
			locP4Constraint->dConstrainMassFlag = false;
			if(locInclusiveP4Constraint == NULL)
			{
				if(dDebugLevel > 5)
					cout << "p4 constraint with pid " << locKinFitParticle->Get_PID() << " as parent marked as full-p4 constraint" << endl;
				locP4Constraint->dIsActualP4ConstraintFlag = true;
				return true; //is not inclusive: we're done
			}
			if(dDebugLevel > 5)
				cout << "p4 constraint with pid " << locKinFitParticle->Get_PID() << " as parent marked as no-mass constraint" << endl;
		}
		else if((locInclusiveP4Constraint != NULL) && (!locUnknownParentFoundFlag)) //if unknown parent, all will be invariant mass constraints
			Mark_AsMissingMassConstraintIfNecessary(locP4Constraint);
	}

	if(locInclusiveP4Constraint == NULL)
		return (!locAreP4ConstraintsFlag); //is not inclusive: if no p4 fit then bad
	return true;
}

void DKinFitter::Mark_AsMissingMassConstraintIfNecessary(DKinFitConstraint_P4* locP4Constraint) const
{
	if(!locP4Constraint->Get_ConstrainMassByInvariantMassFlag())
		return; //already marked as missing mass constraint

	for(size_t loc_i = 0; loc_i < locP4Constraint->dFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locP4Constraint->dFinalParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if(locKinFitParticleType != d_DecayingParticle)
			continue;

		//decaying particle: dive to the next constraint
		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->dP4Constraints;
		for(size_t loc_j = 0; loc_j < locP4Constraints.size(); ++loc_j)
		{
			DKinFitConstraint_P4* locP4SubConstraint = locP4Constraints[loc_j];
			if(locP4SubConstraint == locP4Constraint)
				continue;
			if(!locP4SubConstraint->Get_ConstrainMassByInvariantMassFlag())
			{
				//this next constraint is already by missing mass because X is a decay product of it. therefore, this constraint must be by missing mass also
				locP4Constraint->dConstrainMassByInvariantMassFlag = false;
				if(dDebugLevel > 5)
					cout << "p4 constraint with pid " << locP4Constraint->dInitialParticles[0]->Get_PID() << " as parent marked as a missing-mass constraint" << endl;
			}
			else
			{
				Mark_AsMissingMassConstraintIfNecessary(locP4SubConstraint);
				if(!locP4SubConstraint->Get_ConstrainMassByInvariantMassFlag())
				{
					//this next constraint is already by missing mass because X is a decay product of it. therefore, this constraint must be by missing mass also
					locP4Constraint->dConstrainMassByInvariantMassFlag = false;
					if(dDebugLevel > 5)
						cout << "p4 constraint with pid " << locP4Constraint->dInitialParticles[0]->Get_PID() << " as parent marked as a missing-mass constraint" << endl;
				}
			}
		}
	}
}

bool DKinFitter::Group_Constraints(const deque<DKinFitConstraint_VertexBase*>& locSortedVertexConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints) const
{
	//loop through the vertex constraints
	locSortedConstraints.clear();
	set<DKinFitConstraint_P4*> locHandledP4Constraints;
	for(size_t loc_i = 0; loc_i < locSortedVertexConstraints.size(); ++loc_i)
	{
		//find all p4 constraints that the constraining particles in this vertex constraint are constrained to
		//first get particles
		deque<const DKinFitParticle*> locConstrainParticles = locSortedVertexConstraints[loc_i]->Get_FullConstrainParticles();
		//if spacetime fit, get neutral showers as well
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(locSortedVertexConstraints[loc_i]);
		if(locSpacetimeConstraint != NULL)
		{
			deque<const DKinFitParticle*> locTimeConstrainParticles = locSpacetimeConstraint->Get_OnlyConstrainTimeParticles();
			locConstrainParticles.insert(locConstrainParticles.end(), locTimeConstrainParticles.begin(), locTimeConstrainParticles.end());
		}
		//then get unhandled p4 constraints
		set<DKinFitConstraint_P4*> locP4Constraints;
		for(size_t loc_j = 0; loc_j < locConstrainParticles.size(); ++loc_j)
		{
			DKinFitConstraint_P4* locConstraint = const_cast<DKinFitConstraint_P4*>(locConstrainParticles[loc_j]->Get_ConstrainedAtP4Constraint());
			if(locConstraint != NULL)
				locP4Constraints.insert(locConstraint);
		}
		pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > locConstraintPair(locSortedVertexConstraints[loc_i], locP4Constraints);
		locSortedConstraints.push_back(locConstraintPair);
	}

	return true;
}

bool DKinFitter::Fit_Reaction(void)
{
#ifdef VTRACE
	VT_TRACER("DKinFitter::Fit_Reaction()");
#endif
	if(!Resolve_Constraints())
	{
		dKinFitStatus = d_KinFitFailedSetup;
		return false;
	}

	Set_MatrixSizes();
	Resize_Matrices();
	Fill_InputMatrices();

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
	double locPreviousChiSq = 0.0;
	int locNumIterations = -1;
	TMatrixD locR(dNumF, 1);
	do
	{
		++locNumIterations;
		if(locNumIterations >= int(dMaxNumIterations))
		{
			//sometimes the chisq will walk (very slightly) forever, without any real meaningful change in the variables
			if(dDebugLevel > 10)
				cout << "DKinFitter: At maximum number of iterations, this chisq, last chisq, last resort cutoff = " << dChiSq << ", " << locPreviousChiSq << ", " << dConvergenceChiSqDiff_LastResort << endl;
			if((fabs(dChiSq - locPreviousChiSq) <= dConvergenceChiSqDiff_LastResort) && (dChiSq >= 0.0))
				break; //close enough
			if(dDebugLevel > 10)
				cout << "DKinFitter: Exceeded maximum number of iterations.  Returning false." << endl;
			dKinFitStatus = d_KinFitTooManyIterations;
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
			dKinFitStatus = d_KinFitFailedInversion;
			return false; // matrix is not invertible
		}

		if(dNumXi > 0)
		{
			if(!Calc_dU())
			{
				if(dDebugLevel > 10)
					cout << "DKinFitter: Failed VXi-matrix inversion. Returning false." << endl;
				dKinFitStatus = d_KinFitFailedInversion;
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
	while((fabs(dChiSq - locPreviousChiSq) > dConvergenceChiSqDiff) || (dChiSq < 0.0));

	// dVXi
	if(dNumXi > 0)
		*dVXi = dU;

	// dVEta & dV
	TMatrixDSym locG = dS_Inverse;
	locG.SimilarityT(dF_dEta);
	if(dNumXi > 0)
	{
		TMatrixD locH = dF_dEta_T*dS_Inverse*dF_dXi;
		TMatrixDSym locTempMatrix11 = *dVXi;
		*dVEta = dVY - (locG - locTempMatrix11.Similarity(locH)).Similarity(dVY);

		//dV:
		TMatrixD locEtaXiCovariance = -1.0*dVY*locH*dU;
		for(unsigned int loc_i = 0; loc_i < dNumEta; ++loc_i)
		{
			for(unsigned int loc_j = 0; loc_j < dNumEta; ++loc_j)
				(*dV)(loc_i, loc_j) = (*dVEta)(loc_i, loc_j);
		}
		for(unsigned int loc_i = 0; loc_i < dNumXi; ++loc_i)
		{
			for(unsigned int loc_j = 0; loc_j < dNumXi; ++loc_j)
				(*dV)(loc_i + dNumEta, loc_j + dNumEta) = (*dVXi)(loc_i, loc_j);
		}
		for(unsigned int loc_i = 0; loc_i < dNumEta; ++loc_i)
		{
			for(unsigned int loc_j = 0; loc_j < dNumXi; ++loc_j)
			{
				(*dV)(loc_i, loc_j + dNumEta) = locEtaXiCovariance(loc_i, loc_j);
				(*dV)(loc_j + dNumEta, loc_i) = locEtaXiCovariance(loc_i, loc_j);
			}
		}
	}
	else
	{
		*dVEta = dVY - locG.Similarity(dVY); //destroys locG, but it's not needed anymore
		*dV = *dVEta;
	}

	dEpsilon = dY - dEta;

	Calc_Pulls();
	dNDF = dNumF - dNumXi;
	dConfidenceLevel = TMath::Prob(dChiSq, dNDF);

	Set_FinalTrackInfo();

	if(dDebugLevel > 5)
		cout << "DKinFitter: Final dChiSq, dNDF, dConfidenceLevel = " << dChiSq << ", " << dNDF << ", " << dConfidenceLevel << endl;

	dKinFitStatus = d_KinFitSuccessful;
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
	cout << "DKinFitter: Particle PID, Q, Mass, E, P3, V3, T = " << locKinFitParticle->Get_PID() << ", " << locCharge << ", " << locMass << ", " << locP4.E() << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;
	if(locCovarianceMatrix != NULL)
	{
		cout << "DKinFitter: CovMatrix Diagonal Terms: ";
		for(int loc_i = 0; loc_i < locCovarianceMatrix->GetNcols(); ++loc_i)
			cout << (*locCovarianceMatrix)(loc_i, loc_i) << ", ";
		cout << endl;
	}
	cout << "DKinFitter: Particle E, Px, Vx, Common Vx, T, Common T, L indices = " << locKinFitParticle->Get_EParamIndex() << ", " << locKinFitParticle->Get_PxParamIndex() << ", " << locKinFitParticle->Get_VxParamIndex() << ", " << locKinFitParticle->Get_CommonVxParamIndex() << ", " << locKinFitParticle->Get_TParamIndex() << ", " << locKinFitParticle->Get_CommonTParamIndex() << ", " << locKinFitParticle->Get_LParamIndex() << endl;
	cout << "DKinFitter: Particle CovMatrix E, Px, Vx, T indices = " << locKinFitParticle->Get_CovMatrixEParamIndex() << ", " << locKinFitParticle->Get_CovMatrixPxParamIndex() << ", " << locKinFitParticle->Get_CovMatrixVxParamIndex() << ", " << locKinFitParticle->Get_CovMatrixTParamIndex() << endl;
}

void DKinFitter::Set_MatrixSizes(void)
{
	//set matrix sizes
	dNumXi = 0; //num unknowns
	dNumEta = 0; //num measurables
	dNumF = 0; //num constraint eqs

	//Calculate dNumEta
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = dKinFitParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		if(!locKinFitParticle->Get_IsNeutralShowerFlag())
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) || (locKinFitParticle->Get_ConstrainedAtVertexAndOrTimeConstraint() != NULL))
				dNumEta += 3; //p3
			if(locKinFitParticle->Get_ConstrainedAtVertexAndOrTimeConstraint() != NULL)
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
			if(locKinFitConstraint_P4->Get_IsActualP4ConstraintFlag())
				dNumF += 4; //p4 constraint
			else if(locKinFitConstraint_P4->Get_ConstrainInitialParticleMassFlag())
				dNumF += 1; //mass constraint

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
				if(((locKinFitParticleType == d_MissingParticle) && (locFinalParticles[loc_j]->Get_PID() != 0)) || ((locKinFitParticleType == d_DecayingParticle) && (locFinalParticles[loc_j]->Get_NumP4Constraints() == 1)))
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
			dNumF += 2*locKinFitConstraint_Vertex->dFullConstrainParticles.size();
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->dFullConstrainParticles.size(); ++loc_j)
					cout << locKinFitConstraint_Vertex->dFullConstrainParticles[loc_j]->Get_Charge() << ", " << locKinFitConstraint_Vertex->dFullConstrainParticles[loc_j]->Get_Mass() << "; ";
				cout << endl;
			}
			continue;
		}
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			dNumXi += 4; //v3, t
			deque<DKinFitParticle*> locFullConstrainParticles = locKinFitConstraint_Spacetime->dFullConstrainParticles;
			for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
				dNumF += 3;
			if(dDebugLevel > 10)
			{
				cout << "q's, masses of spacetime vertex constraining particles: ";
				for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
					cout << locFullConstrainParticles[loc_j]->Get_Charge() << ", " << locFullConstrainParticles[loc_j]->Get_Mass() << ";";
				cout << endl;
			}
			if(Get_IsBFieldNearBeamline())
			{
				size_t locNumChargedConstraintParticles = 0;
				size_t locNumDecayingChargedConstraintParticles = 0;
				for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
				{
					if(locFullConstrainParticles[loc_j]->Get_Charge() == 0)
						continue;
					++locNumChargedConstraintParticles;
					if(locFullConstrainParticles[loc_j]->Get_KinFitParticleType() == d_DecayingParticle)
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

	dV->ResizeTo(dNumEta + dNumXi, dNumEta + dNumXi);

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
	dV->Zero();
}

void DKinFitter::Fill_InputMatrices(void)
{
	//fill dY, dEta, dVY, dXi

	DKinFitParticle* locKinFitParticle;
	DKinFitParticleType locKinFitParticleType;
	int locParamIndex, locConstraintIndex_Eta, locConstraintIndex_Xi;
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

		if(!locKinFitParticle->Get_IsNeutralShowerFlag()) //non-neutral shower
		{
			if((locKinFitParticle->Get_IsInP4FitFlag()) || (locKinFitParticle->Get_ConstrainedAtVertexAndOrTimeConstraint() != NULL)) //p3
			{
				locKinFitParticle->Set_PxParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locMomentum.Px();
				dY(locParamIndex + 1, 0) = locMomentum.Py();
				dY(locParamIndex + 2, 0) = locMomentum.Pz();
				locParamIndex += 3;
			}
			if(locKinFitParticle->Get_ConstrainedAtVertexAndOrTimeConstraint() != NULL)
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
	deque<DKinFitParticle*> locNoConstrainParticles;
	deque<DKinFitParticle*> locConstrainParticles;
	bool locRFTimeConstrainedFlag = false;

	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_P4 != NULL)
		{
			locKinFitConstraint_P4->Set_FIndex(locConstraintIndex_Eta);

			if(locKinFitConstraint_P4->Get_IsActualP4ConstraintFlag())
				locConstraintIndex_Eta += 4; //p4
			else if(locKinFitConstraint_P4->Get_ConstrainInitialParticleMassFlag())
				locConstraintIndex_Eta += 1; //mass
			else
				locKinFitConstraint_P4->Set_FIndex(-1);

			DKinFitParticle* locConstrainedKinFitParticle = locKinFitConstraint_P4->dConstrainedP4Particle;
			if(locConstrainedKinFitParticle == NULL)
				continue; //no missing particles or no unknown p3 (p3 is derivable from known p3's (p4 constraint already applied elsewhere))
			DKinFitParticleType locKinFitParticleType = locConstrainedKinFitParticle->Get_KinFitParticleType();
			if(((locKinFitParticleType == d_MissingParticle) && (locConstrainedKinFitParticle->Get_PID() != 0)) || ((locKinFitParticleType == d_DecayingParticle) && (locConstrainedKinFitParticle->Get_NumP4Constraints() == 1)))
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
			locConstrainParticles = locKinFitConstraint_Vertex->dFullConstrainParticles;
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

			locConstrainParticles = locKinFitConstraint_Spacetime->dFullConstrainParticles;
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

			deque<DKinFitParticle*> locFullConstrainParticles = locKinFitConstraint_Spacetime->dFullConstrainParticles;
			for(size_t loc_j = 0; loc_j < locFullConstrainParticles.size(); ++loc_j)
			{
				if(locFullConstrainParticles[loc_j]->Get_Charge() == 0)
					continue;
				locParamIndex = locFullConstrainParticles[loc_j]->Get_LParamIndex();
				if(locParamIndex > 0)
					locFullConstrainParticles[loc_j]->Set_PathLength(dXi(locParamIndex, 0));
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
			//calc energy of missing unknown particle (if present), and set mass appropriately
			if((locKinFitParticleType == d_MissingParticle) && (locKinFitParticle->Get_PID() == 0))
			{
				DKinFitConstraint_P4* locKinFitConstraint_P4 = const_cast<DKinFitConstraint_P4*>(locKinFitParticle->Get_DefinedAtP4Constraint());
				if(locKinFitConstraint_P4 != NULL)
				{
					TLorentzVector locP4;
					Constrain_Particle(locKinFitParticle, locKinFitConstraint_P4, locP4);
					locKinFitParticle->Set_Mass(locP4.M());
				}
			}
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
		//what if a decay product is a missing particle of unknown mass? it's fine: only setting p3 anyway, which is derived
		//locKinFitParticle is now definitely an enclosed decaying particle
		if(locKinFitParticle->Get_DecayingParticleAtProductionVertexFlag())
			locKinFitParticle->Set_Momentum(Calc_DecayingP4(locKinFitConstraint_P4, false).Vect()); //returns p4 at production vertex
		else
			locKinFitParticle->Set_Momentum(Calc_DecayingP4(locKinFitConstraint_P4, true).Vect()); //returns p4 at decay vertex
	}

	if(dRFTimeParamIndex > 0)
		dRFTime = dEta(dRFTimeParamIndex, 0);
}

TLorentzVector DKinFitter::Calc_DecayingP4(DKinFitConstraint_P4* locP4Constraint, bool locDecayMomentumFlag) const
{
	//uses decay products to calculate decaying particle information
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
			deque<DKinFitConstraint_P4*> locP4Constraints = locFinalParticles[loc_i]->dP4Constraints;
			for(size_t loc_j = 0; loc_j < locP4Constraints.size(); ++loc_j)
			{
				if(locP4Constraints[loc_j] == locP4Constraint)
					continue;
				locP4 += Calc_DecayingP4(locP4Constraints[loc_j], true);
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
			int locFIndex = locKinFitConstraint_P4->Get_FIndex();
			if(dDebugLevel > 10)
				cout << "DKinFitter: F index = " << locFIndex << endl;
			if(locFIndex < 0)
				continue; //e.g. a p4 constraint with a decaying particle as an enclosed parent that has no mass constraint

			if(locKinFitConstraint_P4->Get_IsActualP4ConstraintFlag())
			{
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dInitialParticles.size(); ++loc_j)
				{
					locKinFitParticle = (locKinFitConstraint_P4->dInitialParticles)[loc_j];
					Calc_dF_P4(locKinFitConstraint_P4, locKinFitParticle, true, NULL);
				}
				for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
				{
					locKinFitParticle = (locKinFitConstraint_P4->dFinalParticles)[loc_j];
					Calc_dF_P4(locKinFitConstraint_P4, locKinFitParticle, false, NULL);
				}
			}
			else if(locKinFitConstraint_P4->Get_ConstrainInitialParticleMassFlag())
			{
				//mass constraint
				double locTargetedMass = (locKinFitConstraint_P4->dInitialParticles)[0]->Get_Mass();
				if(locKinFitConstraint_P4->Get_ConstrainMassByInvariantMassFlag())
				{
					//invariant mass constraint: loop over decay products
					if(dDebugLevel > 10)
						cout << "invariant mass constraint: loop over decay products" << endl;
					TLorentzVector locDecayingParticleDerivedP4;
					for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
					{
						locKinFitParticle = (locKinFitConstraint_P4->dFinalParticles)[loc_j];
						locDecayingParticleDerivedP4 += Calc_dF_MassP4(locKinFitConstraint_P4, locKinFitParticle, false, true, NULL);
					}
					if(dDebugLevel > 30)
						cout << "Final decaying pxyzE is: " << locDecayingParticleDerivedP4.Px() << ", " << locDecayingParticleDerivedP4.Py() << ", " << locDecayingParticleDerivedP4.Pz() << ", " << locDecayingParticleDerivedP4.E() << endl;
					dF(locFIndex, 0) += locDecayingParticleDerivedP4.M2() - locTargetedMass*locTargetedMass;
					for(size_t loc_j = 0; loc_j < locKinFitConstraint_P4->dFinalParticles.size(); ++loc_j)
					{
						locKinFitParticle = (locKinFitConstraint_P4->dFinalParticles)[loc_j];
						Calc_dF_MassDerivs(locKinFitConstraint_P4, locKinFitParticle, locDecayingParticleDerivedP4, false, true, NULL);
					}
				}
				else //missing particle with unknown mass is a decay product: constrain by missing mass
				{
					//missing mass constraint: use this particle to get the next
					if(dDebugLevel > 10)
						cout << "missing mass constraint: use this particle to get the next" << endl;
					locKinFitParticle = locKinFitConstraint_P4->dInitialParticles[0];
					TLorentzVector locDecayingParticleDerivedP4 = Calc_dF_MassP4(locKinFitConstraint_P4, locKinFitParticle, true, false, NULL);
					if(dDebugLevel > 30)
						cout << "Final decaying pxyzE is: " << locDecayingParticleDerivedP4.Px() << ", " << locDecayingParticleDerivedP4.Py() << ", " << locDecayingParticleDerivedP4.Pz() << ", " << locDecayingParticleDerivedP4.E() << endl;
					Calc_dF_MassDerivs(locKinFitConstraint_P4, locKinFitParticle, locDecayingParticleDerivedP4, true, false, NULL);
					dF(locFIndex, 0) += locDecayingParticleDerivedP4.M2() - locTargetedMass*locTargetedMass;
				}
			}
			continue;
		}

		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(dKinFitConstraints[loc_i]);
		if(locKinFitConstraint_Vertex != NULL)
		{
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Vertex->dFullConstrainParticles.size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Vertex->dFullConstrainParticles)[loc_j];
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
			for(size_t loc_j = 0; loc_j < locKinFitConstraint_Spacetime->dFullConstrainParticles.size(); ++loc_j)
			{
				locKinFitParticle = (locKinFitConstraint_Spacetime->dFullConstrainParticles)[loc_j];
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

void DKinFitter::Calc_dF_P4(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, DKinFitConstraint_P4* locKinFitSubConstraint_P4)
{
	int locFIndex = locKinFitConstraint_P4->Get_FIndex();

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
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locKinFitParticleType != d_DecayingParticle) || (locKinFitParticle->Get_NumP4Constraints() == 1))
	{
		//all but enclosed decaying particle: will instead get p4 from decay products (may still need to propagate it below)
		dF(locFIndex, 0) += locSignMultiplier*locP4.E();
		dF(locFIndex + 1, 0) += locSignMultiplier*locP4.Px();
		dF(locFIndex + 2, 0) += locSignMultiplier*locP4.Py();
		dF(locFIndex + 3, 0) += locSignMultiplier*locP4.Pz();
	}

	if(dDebugLevel > 30)
		cout << "q, mass, sign, pxyzE = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << ", " << locSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(locEnoughVertexFitsFlag && locChargedBFieldFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle))
	{
		//fitting vertex of charged track in magnetic field: momentum changes as function of vertex
		//check decaying particles: don't include propagation factor if p4 is used where it is defined (already propagated)
		//if initial particle is a "detected" particle (actually a decaying particle treated as detected): still propagate vertex (assume p3/v3 defined at production vertex)

		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
		if(dDebugLevel > 30)
			cout << "propagate pxyz by: " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.X() << ", " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.Y() << ", " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.Z() << endl;
		dF(locFIndex + 1, 0) -= locSignMultiplier*locA*locDeltaXCrossH.X();
		dF(locFIndex + 2, 0) -= locSignMultiplier*locA*locDeltaXCrossH.Y();
		dF(locFIndex + 3, 0) -= locSignMultiplier*locA*locDeltaXCrossH.Z();
	}

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

		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->dP4Constraints;
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
				Calc_dF_P4(locKinFitConstraint_P4, locInitialParticles[loc_j], true, locNewKinFitSubConstraint_P4); //else !locInitialStateFlag
			}
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				if(locFinalParticles[loc_j] == locKinFitParticle)
					continue;
				if(dDebugLevel > 30)
					cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
				Calc_dF_P4(locKinFitConstraint_P4, locFinalParticles[loc_j], false, locNewKinFitSubConstraint_P4);
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

TLorentzVector DKinFitter::Calc_dF_MassP4(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, bool locIsInvariantMassConstraint, DKinFitConstraint_P4* locKinFitSubConstraint_P4)
{
	//locIsInvariantMassConstraint: true/false if constraining invariant/missing mass
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
	double locSignMultiplier = (locInitialStateFlag != locIsInvariantMassConstraint) ? 1.0 : -1.0;

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();

	TLorentzVector locP4Sum;
	if((locKinFitParticleType != d_DecayingParticle) || (locKinFitParticle->Get_NumP4Constraints() == 1))
		locP4Sum += locSignMultiplier*locP4; //all but enclosed decaying particle: will instead get p4 from decay products (may still need to propagate it below)

	if(dDebugLevel > 30)
		cout << "q, mass, sign, pxyzE = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << ", " << locSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(locEnoughVertexFitsFlag && locChargedBFieldFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle) && ((locKinFitParticleType != d_DecayingParticle) || (locKinFitSubConstraint_P4 != NULL)))
	{
		//fitting vertex of charged track in magnetic field: momentum changes as function of vertex
		//check decaying particles: don't include propagation factor if p4 is used where it is defined (already propagated)
		//if initial particle is a "detected" particle (actually a decaying particle treated as detected): still propagate vertex (assume p3/v3 defined at production vertex)

		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
		if(dDebugLevel > 30)
			cout << "propagate pxyz by: " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.X() << ", " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.Y() << ", " << -1.0*locSignMultiplier*locA*locDeltaXCrossH.Z() << endl;

		locP4Sum.SetVect(locP4Sum.Vect() - locSignMultiplier*locA*locDeltaXCrossH);
	}

	if((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex < 0))
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassP4() Decaying Particle; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		//get the other p4 constraint its in and loop over the particles in it
		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->dP4Constraints;
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
				locP4Sum += Calc_dF_MassP4(locKinFitConstraint_P4, locInitialParticles[loc_j], true, locIsInvariantMassConstraint, locNewKinFitSubConstraint_P4); //else !locInitialStateFlag
				if(dDebugLevel > 30)
					cout << "pxyzE sum is: " << locP4Sum.Px() << ", " << locP4Sum.Py() << ", " << locP4Sum.Pz() << ", " << locP4Sum.E() << endl;
			}
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				if(locFinalParticles[loc_j] == locKinFitParticle)
					continue;
				if(dDebugLevel > 30)
					cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
				locP4Sum += Calc_dF_MassP4(locKinFitConstraint_P4, locFinalParticles[loc_j], false, locIsInvariantMassConstraint, locNewKinFitSubConstraint_P4);
				if(dDebugLevel > 30)
					cout << "pxyzE sum is: " << locP4Sum.Px() << ", " << locP4Sum.Py() << ", " << locP4Sum.Pz() << ", " << locP4Sum.E() << endl;
			}
		}
	}
	return locP4Sum;
}

void DKinFitter::Calc_dF_MassDerivs(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, TLorentzVector locDecayingParticleDerivedP4, bool locInitialStateFlag, bool locIsInvariantMassConstraint, DKinFitConstraint_P4* locKinFitSubConstraint_P4)
{
	int locFIndex = locKinFitConstraint_P4->Get_FIndex();

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
//	double locSignMultiplier = locInitialStateFlag ? 1.0 : -1.0;
	double locSignMultiplier = (locInitialStateFlag != locIsInvariantMassConstraint) ? 1.0 : -1.0;
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	TVector3 locPXCrossHi = locDecayingParticleDerivedP4.Vect().Cross(locH);

	if(dDebugLevel > 30)
		cout << "q, mass, sign, pxyzE = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << ", " << locSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 1; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		return; //target params are fixed: no partial derivatives
	}
	else if(locChargedBFieldFlag && locEnoughVertexFitsFlag && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle) & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 2; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*locSignMultiplier*(locP4.Px()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*locSignMultiplier*(locP4.Py()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*locSignMultiplier*(locP4.Pz()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Pz());

		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locSignMultiplier*locA*locPXCrossHi.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locSignMultiplier*locA*locPXCrossHi.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locSignMultiplier*locA*locPXCrossHi.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 3; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dF_dEta(locFIndex, locEParamIndex) = 2.0*locSignMultiplier*(locDecayingParticleDerivedP4.E() - locEOverPSq*locDecayingParticleDerivedP4.Vect().Dot(locP4.Vect()));

		double locDeltaXDotPXOverMagDeltaXSq = locDeltaX.Dot(locDecayingParticleDerivedP4.Vect())/(locDeltaX.Mag2());
		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locSignMultiplier*locP4.Px()*(locDecayingParticleDerivedP4.Px()/locDeltaX.X() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locSignMultiplier*locP4.Py()*(locDecayingParticleDerivedP4.Py()/locDeltaX.Y() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locSignMultiplier*locP4.Pz()*(locDecayingParticleDerivedP4.Pz()/locDeltaX.Z() - locDeltaXDotPXOverMagDeltaXSq);

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex >= 0)))
	{
		//missing or open-ended-decaying particle: p3 is unknown (not derivable) //must be the constrained particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 4; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dXi(locFIndex, locPxParamIndex) = 2.0*locSignMultiplier*(locP4.Px()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Px());
		dF_dXi(locFIndex, locPxParamIndex + 1) = 2.0*locSignMultiplier*(locP4.Py()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Py());
		dF_dXi(locFIndex, locPxParamIndex + 2) = 2.0*locSignMultiplier*(locP4.Pz()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Pz());
	}
	else if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 5; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;
		if((locCharge != 0) && Get_IsBFieldNearBeamline() && (locKinFitParticle->Get_NumVertexFits() == 2) && (locKinFitSubConstraint_P4 != NULL))
		{
			//if locKinFitSubConstraint_P4 is NULL is first particle; don't include: replace it
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_MassDerivs() Section 5a" << endl;

			dF_dXi(locFIndex, locVxParamIndex) += 2.0*locSignMultiplier*locA*locPXCrossHi.X();
			dF_dXi(locFIndex, locVxParamIndex + 1) += 2.0*locSignMultiplier*locA*locPXCrossHi.Y();
			dF_dXi(locFIndex, locVxParamIndex + 2) += 2.0*locSignMultiplier*locA*locPXCrossHi.Z();

			dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);
		}

		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->dP4Constraints;
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
				Calc_dF_MassDerivs(locKinFitConstraint_P4, locInitialParticles[loc_j], locDecayingParticleDerivedP4, true, locIsInvariantMassConstraint, locNewKinFitSubConstraint_P4); //else !locInitialStateFlag
			}
			for(size_t loc_j = 0; loc_j < locFinalParticles.size(); ++loc_j)
			{
				if(locFinalParticles[loc_j] == locKinFitParticle)
					continue;
				if(dDebugLevel > 30)
					cout << "decaying, partially replace with final-state q, mass = " << locFinalParticles[loc_j]->Get_Charge() << ", " << locFinalParticles[loc_j]->Get_Mass() << endl;
				Calc_dF_MassDerivs(locKinFitConstraint_P4, locFinalParticles[loc_j], locDecayingParticleDerivedP4, false, locIsInvariantMassConstraint, locNewKinFitSubConstraint_P4);
			}
		}
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 6; q, mass = " << locKinFitParticle->Get_Charge() << ", " << locKinFitParticle->Get_Mass() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*locSignMultiplier*(locP4.Px()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*locSignMultiplier*(locP4.Py()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*locSignMultiplier*(locP4.Pz()*locDecayingParticleDerivedP4.E()/locP4.E() - locDecayingParticleDerivedP4.Pz());
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
		//replace the momentum with it's component tracks
		deque<DKinFitConstraint_P4*> locP4Constraints = locKinFitParticle->dP4Constraints;
		DKinFitConstraint_P4* locKinFitConstraint_P4 = NULL;

		//true if the object's p3, v3, & t are defined at its production vertex (when its a final state particle). else at it's decay vertex (when initial state)
		bool locAtProductionVertexFlag = locKinFitParticle->Get_DecayingParticleAtProductionVertexFlag();

		// find the p4 constraint where the decaying particle momentum is defined //this is the point where the position is defined as well
			//why this constraint? because if you choose the other one, the momentum has to be propagated across a distance, requiring a different (& MUCH more complicated) partial derivative
		for(size_t loc_i = 0; loc_i < locP4Constraints.size(); ++loc_i)
		{
			deque<const DKinFitParticle*> locParticlesToSearch = locAtProductionVertexFlag ? locP4Constraints[loc_i]->Get_FinalParticles() : locP4Constraints[loc_i]->Get_InitialParticles();
			bool locMatchFoundFlag = false;
			for(size_t loc_j = 0; loc_j < locParticlesToSearch.size(); ++loc_j)
			{
				if(locParticlesToSearch[loc_j] != locKinFitParticle)
					continue;
				locMatchFoundFlag = true;
				break;
			}
			if(!locMatchFoundFlag)
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

		//decaying particles with momentum derived from other particles
		if((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex < 0))
		{
			//enclosed decaying particle: the momentum is derived from the momentum of its decay products
			deque<DKinFitConstraint_P4*> locConstraints = dKinFitParticles[loc_i]->dP4Constraints;
			DKinFitConstraint_P4* locConstraintAsParent = NULL;
			for(size_t loc_j = 0; loc_j < locConstraints.size(); ++loc_j)
			{
				if((locConstraints[loc_j]->Get_InitialParticles())[0] != dKinFitParticles[loc_i])
					continue;
				locConstraintAsParent = locConstraints[loc_j];
				break;
			}
			if(locConstraintAsParent != NULL)
			{
				//decaying particle used in a p4 constraint
				TMatrixD locJacobian(7, dNumEta + dNumXi);
				locJacobian.Zero();
				if(dKinFitParticles[loc_i]->Get_DecayingParticleAtProductionVertexFlag())
					Calc_DecayingParticleJacobian(locConstraintAsParent, false, locJacobian); //computes jacobian at production vertex
				else
					Calc_DecayingParticleJacobian(locConstraintAsParent, true, locJacobian); //computes jacobian at decay vertex

				//set vertex & time terms
				if(locVxParamIndex >= 0)
				{
					locJacobian(locCovMatrixVxParamIndex, locVxParamIndex + dNumEta) = 1.0;
					locJacobian(locCovMatrixVxParamIndex + 1, locVxParamIndex + dNumEta + 1) = 1.0;
					locJacobian(locCovMatrixVxParamIndex + 2, locVxParamIndex + dNumEta + 2) = 1.0;
				}
				int locTParamIndex = dKinFitParticles[loc_i]->Get_TParamIndex();
				if(locTParamIndex >= 0)
					locJacobian(6, locTParamIndex + dNumEta) = 1.0;

				TMatrixDSym locTempMatrix = *dV;
				locCovarianceMatrix = locTempMatrix.Similarity(locJacobian);
			}
			else
			{
				//decaying particle not used in a p4 constraint (e.g. vertex-only fit)
				if(locVxParamIndex >= 0)
				{
					for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
					{
						for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
							locCovarianceMatrix(loc_j + locCovMatrixVxParamIndex, loc_k + locCovMatrixVxParamIndex) = (*dVXi)(locVxParamIndex + loc_j, locVxParamIndex + loc_k);
					}
				}
				if(locTParamIndex >= 0)
					locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixTParamIndex) = (*dVXi)(locTParamIndex, locTParamIndex);
				if((locVxParamIndex >= 0) && (locTParamIndex >= 0)) //both included in the fit
				{
					for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
					{
						locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) = (*dVXi)(locVxParamIndex + loc_j, locTParamIndex + 0);
						locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) = (*dVXi)(locTParamIndex + 0, locVxParamIndex + loc_j);
					}
				}
			}

			if(dDebugLevel >= 50)
			{
				cout << "FINAL COV MATRIX (enclosed decaying particle):" << endl;
				Print_Matrix(locCovarianceMatrix);
			}
			continue;
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

		//cross terms: E & V (neutral shower)
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
			double locDenominator = sqrt(dVY(locEParamIndex, locEParamIndex));
			locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/locDenominator : 0.0;
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
				double locDenominator = sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
				locCovarianceMatrix(locCovMatrixEParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixEParamIndex + 0) *= locUncertaintyRatio;
			}
		}

		//cross terms: E & T (neutral shower)
		if((locEParamIndex >= 0) && (locTParamIndex >= 0)) //both included in the fit
		{
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) = locKinFitMatrix(locEParamIndex, locTParamIndex);
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) = locKinFitMatrix(locTParamIndex, locEParamIndex);
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex >= 0) && (locTParamIndex < 0)) //only E included in the fit
		{
			double locDenominator = sqrt(dVY(locEParamIndex, locEParamIndex));
			locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/locDenominator : 0.0;
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) *= locUncertaintyRatio;
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) *= locUncertaintyRatio;
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex < 0) && (locTParamIndex >= 0) && (locCovMatrixEParamIndex >= 0)) //only T included in the fit //E may not be in the covariance matrix!!
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) *= locUncertaintyRatio;
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) *= locUncertaintyRatio;
		}

		//cross terms: P & V
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
				double locDenominator = sqrt(dVY(locPxParamIndex + loc_j, locPxParamIndex + loc_j));
				locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/locDenominator : 0.0;
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
				double locDenominator = sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
				for(unsigned int loc_k = 0; loc_k < 3; ++loc_k)
				{
					locCovarianceMatrix(locCovMatrixPxParamIndex + loc_k, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
					locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixPxParamIndex + loc_k) *= locUncertaintyRatio;
				}
			}
		}

		//cross terms: P & T
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
				double locDenominator = sqrt(dVY(locPxParamIndex + loc_j, locPxParamIndex + loc_j));
				locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/locDenominator : 0.0;
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}

		//cross terms: V & T
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
				double locDenominator = sqrt(dVY(locVxParamIndex + loc_j, locVxParamIndex + loc_j));
				locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locVxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
			for(unsigned int loc_j = 0; loc_j < 3; ++loc_j)
			{
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}

		if(dDebugLevel >= 50)
		{
			cout << "FINAL COV MATRIX:" << endl;
			Print_Matrix(locCovarianceMatrix);
		}
	} //end set cov matrix loop

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

void DKinFitter::Calc_DecayingParticleJacobian(DKinFitConstraint_P4* locP4Constraint, bool locDecayVertexFlag, TMatrixD& locJacobian) const
{
	//locJacobian: matrix used to convert dV to the decaying particle covariance matrix: indices are px, py, pz, x, y, z, t
		//dimensions are: 7, (dNumXi + dNumEta);
	//uses decay products to calculate decaying particle information
	//locDecayVertexFlag = true to compute jacobian at the decay vertex, false at the production vertex

	if(!locDecayVertexFlag)
	{
		if(dDebugLevel > 50)
			cout << "compute jacobian at production vertex" << endl;
		//propagate decaying particle momentum from the decay vertex to the production vertex (if necessary: charged, b-field, etc.)
		deque<DKinFitParticle*> locInitialParticles = locP4Constraint->dInitialParticles;
		DKinFitParticle* locKinFitParticle = locInitialParticles[0]; //the decaying particle

		size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
		int locCharge = locKinFitParticle->Get_Charge();
		bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();

		if(locChargedBFieldFlag && (locNumVertexFits == 2))
		{
			if(dDebugLevel > 50)
				cout << "charged, enclosed decaying particle in a b-field in vertex fits" << endl;
			//charged, enclosed decaying particle in a b-field in vertex fits
			TVector3 locPosition = locKinFitParticle->Get_Position();
			TVector3 locBField = Get_BField(locPosition);
			TVector3 locH = locBField.Unit();
			double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

			int locVxParamIndex = locKinFitParticle->Get_VxParamIndex() + dNumEta;
			int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex() + dNumEta;

			//vertex factors
			locJacobian(0, locVxParamIndex + 1) += locA*locH.Z();
			locJacobian(0, locVxParamIndex + 2) += -1.0*locA*locH.Y();

			locJacobian(1, locVxParamIndex) += -1.0*locA*locH.Z();
			locJacobian(1, locVxParamIndex + 2) += locA*locH.X();

			locJacobian(2, locVxParamIndex) += locA*locH.Y();
			locJacobian(2, locVxParamIndex + 1) += -1.0*locA*locH.X();

			locJacobian(0, locCommonVxParamIndex + 1) -= locJacobian(0, locVxParamIndex + 1);
			locJacobian(0, locCommonVxParamIndex + 2) -= locJacobian(0, locVxParamIndex + 2);

			locJacobian(1, locCommonVxParamIndex) -= locJacobian(1, locVxParamIndex);
			locJacobian(1, locCommonVxParamIndex + 2) -= locJacobian(1, locVxParamIndex + 2);

			locJacobian(2, locCommonVxParamIndex) -= locJacobian(2, locVxParamIndex);
			locJacobian(2, locCommonVxParamIndex + 1) -= locJacobian(2, locVxParamIndex + 1);
		}
	}

	deque<DKinFitParticle*> locFinalParticles = locP4Constraint->dFinalParticles;
	for(size_t loc_i = 0; loc_i < locFinalParticles.size(); ++loc_i)
	{
		DKinFitParticle* locKinFitParticle = locFinalParticles[loc_i];
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		int locCharge = locKinFitParticle->Get_Charge();
		TLorentzVector locP4 = locKinFitParticle->Get_P4();
		TVector3 locPosition = locKinFitParticle->Get_Position();
		TVector3 locBField = Get_BField(locPosition);
		TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
		TVector3 locDeltaX = locCommonVertex - locPosition;
		if(dDebugLevel > 50)
			cout << "jacobian: decay product: q, mass = " << locCharge << ", " << locKinFitParticle->Get_Mass() << endl;

		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		size_t locNumVertexFits = locKinFitParticle->Get_NumVertexFits();
		bool locEnoughVertexFitsFlag = (locNumVertexFits > 0) && ((locNumVertexFits == 2) || (locKinFitParticleType != d_DecayingParticle));
		bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();
		bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

		int locEParamIndex = locKinFitParticle->Get_EParamIndex();
		int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locPxParamIndex += dNumEta;
		int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			locVxParamIndex += dNumEta;
		int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex() + dNumEta;

		if(locKinFitParticleType == d_TargetParticle)
			continue;
		else if(locChargedBFieldFlag && locEnoughVertexFitsFlag && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
		{
			if(dDebugLevel > 50)
				cout << "jacobian: partials part 1" << endl;

			locJacobian(0, locPxParamIndex) = 1.0;
			locJacobian(1, locPxParamIndex + 1) = 1.0;
			locJacobian(2, locPxParamIndex + 2) = 1.0;

			locJacobian(0, locVxParamIndex + 1) = locA*locH.Z();
			locJacobian(0, locVxParamIndex + 2) = -1.0*locA*locH.Y();

			locJacobian(1, locVxParamIndex) = -1.0*locA*locH.Z();
			locJacobian(1, locVxParamIndex + 2) = locA*locH.X();

			locJacobian(2, locVxParamIndex) = locA*locH.Y();
			locJacobian(2, locVxParamIndex + 1) = -1.0*locA*locH.X();

			locJacobian(0, locCommonVxParamIndex + 1) -= locJacobian(0, locVxParamIndex + 1);
			locJacobian(0, locCommonVxParamIndex + 2) -= locJacobian(0, locVxParamIndex + 2);

			locJacobian(1, locCommonVxParamIndex) -= locJacobian(1, locVxParamIndex);
			locJacobian(1, locCommonVxParamIndex + 2) -= locJacobian(1, locVxParamIndex + 2);

			locJacobian(2, locCommonVxParamIndex) -= locJacobian(2, locVxParamIndex);
			locJacobian(2, locCommonVxParamIndex + 1) -= locJacobian(2, locVxParamIndex + 1);
		}
		else if(locNeutralShowerFlag)
		{
			if(dDebugLevel > 50)
				cout << "jacobian: partials part 2" << endl;

			double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
			locJacobian(0, locEParamIndex) = locEOverPSq*locP4.Px();
			locJacobian(1, locEParamIndex) = locEOverPSq*locP4.Py();
			locJacobian(2, locEParamIndex) = locEOverPSq*locP4.Pz();

			TVector3 locDeltaXOverMagDeltaXSq = locDeltaX*(1.0/locDeltaX.Mag2());

			locJacobian(0, locVxParamIndex) = locP4.Px()*(locDeltaXOverMagDeltaXSq.X() - 1.0/locDeltaX.X());
			locJacobian(1, locVxParamIndex + 1) = locP4.Py()*(locDeltaXOverMagDeltaXSq.Y() - 1.0/locDeltaX.Y());
			locJacobian(2, locVxParamIndex + 2) = locP4.Pz()*(locDeltaXOverMagDeltaXSq.Z() - 1.0/locDeltaX.Z());

			locJacobian(0, locVxParamIndex + 1) = locP4.Px()*locDeltaXOverMagDeltaXSq.Y();
			locJacobian(0, locVxParamIndex + 2) = locP4.Px()*locDeltaXOverMagDeltaXSq.Z();

			locJacobian(1, locVxParamIndex) = locP4.Py()*locDeltaXOverMagDeltaXSq.X();
			locJacobian(1, locVxParamIndex + 2) = locP4.Py()*locDeltaXOverMagDeltaXSq.Z();

			locJacobian(2, locVxParamIndex) = locP4.Pz()*locDeltaXOverMagDeltaXSq.X();
			locJacobian(2, locVxParamIndex + 1) = locP4.Pz()*locDeltaXOverMagDeltaXSq.Y();

			locJacobian(0, locCommonVxParamIndex) -= locJacobian(0, locVxParamIndex);
			locJacobian(1, locCommonVxParamIndex + 1) -= locJacobian(1, locVxParamIndex + 1);
			locJacobian(2, locCommonVxParamIndex + 2) -= locJacobian(2, locVxParamIndex + 2);

			locJacobian(0, locCommonVxParamIndex + 1) -= locJacobian(0, locVxParamIndex + 1);
			locJacobian(0, locCommonVxParamIndex + 2) -= locJacobian(0, locVxParamIndex + 2);

			locJacobian(1, locCommonVxParamIndex) -= locJacobian(1, locVxParamIndex);
			locJacobian(1, locCommonVxParamIndex + 2) -= locJacobian(1, locVxParamIndex + 2);

			locJacobian(2, locCommonVxParamIndex) -= locJacobian(2, locVxParamIndex);
			locJacobian(2, locCommonVxParamIndex + 1) -= locJacobian(2, locVxParamIndex + 1);
		}
		else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex >= 0)))
		{
			if(dDebugLevel > 50)
				cout << "jacobian: partials part 3" << endl;

			//missing or open-ended-decaying particle: p3 is unknown (not derivable)
			locJacobian(0, locPxParamIndex) = 1.0;
			locJacobian(1, locPxParamIndex + 1) = 1.0;
			locJacobian(2, locPxParamIndex + 2) = 1.0;
		}
		else if(locKinFitParticleType == d_DecayingParticle)
		{
			if(dDebugLevel > 50)
				cout << "jacobian: partials part 4" << endl;

			//charged, enclosed decaying particle in a b-field
			if(locChargedBFieldFlag && (locKinFitParticle->Get_NumVertexFits() == 2))
			{
				if(dDebugLevel > 50)
					cout << "jacobian: partials part 4a" << endl;

				//vertex factors
				locJacobian(0, locVxParamIndex + 1) += locA*locH.Z();
				locJacobian(0, locVxParamIndex + 2) += -1.0*locA*locH.Y();

				locJacobian(1, locVxParamIndex) += -1.0*locA*locH.Z();
				locJacobian(1, locVxParamIndex + 2) += locA*locH.X();

				locJacobian(2, locVxParamIndex) += locA*locH.Y();
				locJacobian(2, locVxParamIndex + 1) += -1.0*locA*locH.X();

				locJacobian(0, locCommonVxParamIndex + 1) -= locJacobian(0, locVxParamIndex + 1);
				locJacobian(0, locCommonVxParamIndex + 2) -= locJacobian(0, locVxParamIndex + 2);

				locJacobian(1, locCommonVxParamIndex) -= locJacobian(1, locVxParamIndex);
				locJacobian(1, locCommonVxParamIndex + 2) -= locJacobian(1, locVxParamIndex + 2);

				locJacobian(2, locCommonVxParamIndex) -= locJacobian(2, locVxParamIndex);
				locJacobian(2, locCommonVxParamIndex + 1) -= locJacobian(2, locVxParamIndex + 1);
			}

			//p4 is derived from other particles:
			deque<DKinFitConstraint_P4*> locP4Constraints = locFinalParticles[loc_i]->dP4Constraints;
			for(size_t loc_j = 0; loc_j < locP4Constraints.size(); ++loc_j)
			{
				if(locP4Constraints[loc_j] == locP4Constraint)
					continue;
				if(dDebugLevel > 50)
					cout << "jacobian: partials part 4b" << endl;
				Calc_DecayingParticleJacobian(locP4Constraints[loc_j], true, locJacobian);
				break;
			}
		}
		else
		{
			if(dDebugLevel > 50)
				cout << "jacobian: partials part 5" << endl;

			// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
			locJacobian(0, locPxParamIndex) = 1.0;
			locJacobian(1, locPxParamIndex + 1) = 1.0;
			locJacobian(2, locPxParamIndex + 2) = 1.0;
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

	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();
	int locCommonTParamIndex = locKinFitParticle->Get_CommonTParamIndex();

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
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

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

