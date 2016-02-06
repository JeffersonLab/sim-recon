#include "DKinFitUtils.h"

/*************************************************************** RESOURCE MANAGEMENT ***************************************************************/

DKinFitUtils::DKinFitUtils(void)
{
	dKinFitter = NULL; //Is set by DKinFitter constructor
	dLinkVerticesFlag = true;
	dDebugLevel = 0;

	size_t locNumFitsPerEventGuess = 20;

	//final x2: input/output
	dMaxKinFitParticlePoolSize = 5*locNumFitsPerEventGuess*2;

	dMaxKinFitConstraintVertexPoolSize = 2*locNumFitsPerEventGuess*2;
	dMaxKinFitConstraintSpacetimePoolSize = 2*locNumFitsPerEventGuess*2;
	dMaxKinFitConstraintP4PoolSize = locNumFitsPerEventGuess*2;
	dMaxKinFitConstraintMassPoolSize = 2*locNumFitsPerEventGuess*2;

	dMaxKinFitChainPoolSize = locNumFitsPerEventGuess;
	dMaxKinFitChainStepPoolSize = 3*locNumFitsPerEventGuess;

	dMaxMatrixDSymPoolSize = 5*locNumFitsPerEventGuess*2;
	dMaxLargeMatrixDSymPoolSize = 3*locNumFitsPerEventGuess;
}

DKinFitUtils::~DKinFitUtils(void)
{
	for(size_t loc_i = 0; loc_i < dKinFitParticlePool_All.size(); ++loc_i)
		delete dKinFitParticlePool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintVertexPool_All.size(); ++loc_i)
		delete dKinFitConstraintVertexPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintSpacetimePool_All.size(); ++loc_i)
		delete dKinFitConstraintSpacetimePool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintP4Pool_All.size(); ++loc_i)
		delete dKinFitConstraintP4Pool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitConstraintMassPool_All.size(); ++loc_i)
		delete dKinFitConstraintMassPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitChainPool_All.size(); ++loc_i)
		delete dKinFitChainPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dKinFitChainStepPool_All.size(); ++loc_i)
		delete dKinFitChainStepPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dMatrixDSymPool_All.size(); ++loc_i)
		delete dMatrixDSymPool_All[loc_i];

	for(size_t loc_i = 0; loc_i < dLargeMatrixDSymPool_All.size(); ++loc_i)
		delete dLargeMatrixDSymPool_All[loc_i];
}

void DKinFitUtils::Reset_NewEvent(void)
{
	dParticleMap_OutputToInput.clear();
	dMassConstraintMap.clear();
	dP4ConstraintMap.clear();
	dVertexConstraintMap.clear();
	dSpacetimeConstraintMap.clear();

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

	if(dKinFitConstraintMassPool_All.size() > dMaxKinFitConstraintMassPoolSize)
	{
		for(size_t loc_i = dMaxKinFitConstraintMassPoolSize; loc_i < dKinFitConstraintMassPool_All.size(); ++loc_i)
			delete dKinFitConstraintMassPool_All[loc_i];
		dKinFitConstraintMassPool_All.resize(dMaxKinFitConstraintMassPoolSize);
	}
	dKinFitConstraintMassPool_Available = dKinFitConstraintMassPool_All;

	if(dKinFitChainPool_All.size() > dMaxKinFitChainPoolSize)
	{
		for(size_t loc_i = dMaxKinFitChainPoolSize; loc_i < dKinFitChainPool_All.size(); ++loc_i)
			delete dKinFitChainPool_All[loc_i];
		dKinFitChainPool_All.resize(dMaxKinFitChainPoolSize);
	}
	dKinFitChainPool_Available = dKinFitChainPool_All;

	if(dKinFitChainStepPool_All.size() > dMaxKinFitChainStepPoolSize)
	{
		for(size_t loc_i = dMaxKinFitChainStepPoolSize; loc_i < dKinFitChainStepPool_All.size(); ++loc_i)
			delete dKinFitChainStepPool_All[loc_i];
		dKinFitChainStepPool_All.resize(dMaxKinFitChainStepPoolSize);
	}
	dKinFitChainStepPool_Available = dKinFitChainStepPool_All;

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

	Reset_NewFit();
}

void DKinFitUtils::Preallocate_MatrixMemory(void)
{
	//pre-allocate matrix memory
	if(dMatrixDSymPool_Available.empty())
	{
		for(size_t loc_i = 0; loc_i < dMaxMatrixDSymPoolSize; ++loc_i)
		{
			TMatrixDSym* locMatrix = Get_MatrixDSymResource();	
			locMatrix->ResizeTo(7, 7);
		}
		dMatrixDSymPool_Available = dMatrixDSymPool_All;
	}

	//pre-allocate large matrix memory
	if(dLargeMatrixDSymPool_Available.empty())
	{
		for(size_t loc_i = 0; loc_i < dMaxLargeMatrixDSymPoolSize; ++loc_i)
		{
			TMatrixDSym* locMatrix = Get_LargeMatrixDSymResource();	
			locMatrix->ResizeTo(100, 100);
		}
		dLargeMatrixDSymPool_Available = dLargeMatrixDSymPool_All;
	}
}

DKinFitParticle* DKinFitUtils::Get_KinFitParticleResource(void)
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

DKinFitConstraint_Vertex* DKinFitUtils::Get_KinFitConstraintVertexResource(void)
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

DKinFitConstraint_Spacetime* DKinFitUtils::Get_KinFitConstraintSpacetimeResource(void)
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

DKinFitConstraint_P4* DKinFitUtils::Get_KinFitConstraintP4Resource(void)
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

DKinFitConstraint_Mass* DKinFitUtils::Get_KinFitConstraintMassResource(void)
{
	DKinFitConstraint_Mass* locKinFitConstraint;
	if(dKinFitConstraintMassPool_Available.empty())
	{
		locKinFitConstraint = new DKinFitConstraint_Mass;
		dKinFitConstraintMassPool_All.push_back(locKinFitConstraint);
	}
	else
	{
		locKinFitConstraint = dKinFitConstraintMassPool_Available.back();
		locKinFitConstraint->Reset();
		dKinFitConstraintMassPool_Available.pop_back();
	}
	return locKinFitConstraint;
}

DKinFitChain* DKinFitUtils::Get_KinFitChainResource(void)
{
	DKinFitChain* locKinFitChain;
	if(dKinFitChainPool_Available.empty())
	{
		locKinFitChain = new DKinFitChain;
		dKinFitChainPool_All.push_back(locKinFitChain);
	}
	else
	{
		locKinFitChain = dKinFitChainPool_Available.back();
		locKinFitChain->Reset();
		dKinFitChainPool_Available.pop_back();
	}
	return locKinFitChain;
}

DKinFitChainStep* DKinFitUtils::Get_KinFitChainStepResource(void)
{
	DKinFitChainStep* locKinFitChainStep;
	if(dKinFitChainStepPool_Available.empty())
	{
		locKinFitChainStep = new DKinFitChainStep;
		dKinFitChainStepPool_All.push_back(locKinFitChainStep);
	}
	else
	{
		locKinFitChainStep = dKinFitChainStepPool_Available.back();
		locKinFitChainStep->Reset();
		dKinFitChainStepPool_Available.pop_back();
	}
	return locKinFitChainStep;
}

TMatrixDSym* DKinFitUtils::Get_MatrixDSymResource(void)
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

TMatrixDSym* DKinFitUtils::Get_LargeMatrixDSymResource(void)
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

/***************************************************************** CREATE PARTICLES ****************************************************************/

DKinFitParticle* DKinFitUtils::Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex);

	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_BeamParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Beam particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils::Make_TargetParticle(int locPID, int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_TargetParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Target particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils::Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Detected particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils::Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 5) || (locCovarianceMatrix->GetNcols() != 5))
		return NULL; //is not 5x5

	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(0);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_IsNeutralShowerFlag(true);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex); //will be updated with vertex guess
	locKinFitParticle->Set_ShowerEnergy(locShowerEnergy);
	locKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Detected shower created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils::Make_MissingParticle(int locPID, int locCharge, double locMass)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_MissingParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Missing particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

DKinFitParticle* DKinFitUtils::Make_DecayingParticle(int locPID, int locCharge, double locMass, const set<DKinFitParticle*>& locFromInitialState, const set<DKinFitParticle*>& locFromFinalState)
{
	DKinFitParticle* locKinFitParticle = Get_KinFitParticleResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_KinFitParticleType(d_DecayingParticle);
	locKinFitParticle->Set_FromInitialState(locFromInitialState);
	locKinFitParticle->Set_FromFinalState(locFromFinalState);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Decaying particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

/************************************************************* SETUP VERTEX CONSTRAINTS ************************************************************/

deque<DKinFitConstraint_Vertex*> DKinFitUtils::Create_VertexConstraints(const DKinFitChain* locKinFitChain, bool locSpacetimeFitFlag)
{
	deque<set<DKinFitParticle*> > locAllVertices = Setup_VertexConstraints(locKinFitChain);
	return Create_VertexConstraints(locAllVertices, locSpacetimeFitFlag);
}

deque<set<DKinFitParticle*> > DKinFitUtils::Setup_VertexConstraints(const DKinFitChain* locKinFitChain)
{
	//This function organizes the particles into different vertices
	//It does NOT determine whether there is enough information to actually FIT the vertices in question
	if(dDebugLevel > 10)
		cout << "DKinFitUtils: Organize Vertex & Spacetime Constraints" << endl;

	deque<set<DKinFitParticle*> > locAllVertices;
	set<size_t> locIncludedStepIndices;
	for(size_t loc_i = 0; loc_i < locKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
	{
		if(locIncludedStepIndices.find(loc_i) != locIncludedStepIndices.end())
			continue; //already did this step

		//Start a new vertex and save it
		set<DKinFitParticle*> locVertexParticles;
		Setup_VertexConstraint(locKinFitChain, loc_i, locVertexParticles, locIncludedStepIndices);
		locAllVertices.push_back(locVertexParticles);
	}

	if(dDebugLevel > 10)
		cout << "DKinFitUtils: Vertex & spacetime constraints organized" << endl;

	return locAllVertices;
}

void DKinFitUtils::Setup_VertexConstraint(const DKinFitChain* locKinFitChain, size_t locStepIndex, set<DKinFitParticle*>& locVertexParticles, set<size_t>& locIncludedStepIndices)
{
	locIncludedStepIndices.insert(locStepIndex);

	//get the step
	const DKinFitChainStep* locKinFitChainStep = locKinFitChain->Get_KinFitChainStep(locStepIndex);

	//add all particles to the vertex constraint
	set<DKinFitParticle*> locAllParticles = locKinFitChainStep->Get_AllParticles();
	locVertexParticles.insert(locAllParticles.begin(), locAllParticles.end());

	//loop over final particles: dive through decaying particles that decay in-place
	set<DKinFitParticle*> locFinalParticles = locKinFitChainStep->Get_FinalParticles();
	set<DKinFitParticle*>::const_iterator locParticleIterator = locFinalParticles.begin();
	for(; locParticleIterator != locFinalParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		//check if particle decays in-place
		if((locKinFitParticleType == d_DecayingParticle) && !Get_IsDetachedVertex(locKinFitParticle->Get_PID()))
		{
			//yes: combine with the decay products
			int locDecayStepIndex = locKinFitChain->Get_DecayStepIndex(locKinFitParticle);
			Setup_VertexConstraint(locKinFitChain, locDecayStepIndex, locVertexParticles, locIncludedStepIndices);
		}
	}
}

deque<DKinFitConstraint_Vertex*> DKinFitUtils::Create_VertexConstraints(const deque<set<DKinFitParticle*> >& locAllVertices, bool locSpacetimeFitFlag)
{
	if(dDebugLevel > 10)
		cout << "DKinFitUtils: Create vertex constraints." << endl;

	//resolve links between vertex & time fits (decaying particles), create constraints, sort them, and return them
		//sort: the order in which they are defined (as required by the decaying particles they're using)
	deque<DKinFitConstraint_Vertex*> locSortedConstraints;

	//initialize particle groupings
	deque<set<DKinFitParticle*> > locAllFullConstrainParticles, locAllDecayingParticles, locAllOnlyConstrainTimeParticles, locAllNoConstrainParticles;
	for(size_t loc_i = 0; loc_i < locAllVertices.size(); ++loc_i)
	{
		set<DKinFitParticle*> locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles;
		Group_VertexParticles(locAllVertices[loc_i], locFullConstrainParticles, locDecayingParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles);
		locAllFullConstrainParticles.push_back(locFullConstrainParticles);
		locAllDecayingParticles.push_back(locDecayingParticles);
		locAllOnlyConstrainTimeParticles.push_back(locOnlyConstrainTimeParticles);
		locAllNoConstrainParticles.push_back(locNoConstrainParticles);
	}

	//initial pass through vertices:
		//only accept vertices with at least 2 detected, constraining particles
		//will accept decaying particles on the next pass
	set<DKinFitParticle*> locDefinedDecayingParticles;
	for(size_t loc_i = 0; loc_i < locAllFullConstrainParticles.size();)
	{
		if(locAllFullConstrainParticles[loc_i].size() < 2)
		{
			++loc_i;
			continue;
		}

		//Create vertex/spacetime constraint, with decaying particles as no-constrain particles
		locAllNoConstrainParticles[loc_i].insert(locAllDecayingParticles[loc_i].begin(), locAllDecayingParticles[loc_i].end());
		if(locSpacetimeFitFlag)
			locSortedConstraints.push_back(Make_SpacetimeConstraint(locAllFullConstrainParticles[loc_i], locAllOnlyConstrainTimeParticles[loc_i], locAllNoConstrainParticles[loc_i]));
		else //vertex only
		{
			locAllNoConstrainParticles[loc_i].insert(locAllOnlyConstrainTimeParticles[loc_i].begin(), locAllOnlyConstrainTimeParticles[loc_i].end());
			locSortedConstraints.push_back(Make_VertexConstraint(locAllFullConstrainParticles[loc_i], locAllNoConstrainParticles[loc_i]));
		}

		//The positions of these decaying particles are now defined: Can use to constrain vertices in later constraints
		locDefinedDecayingParticles.insert(locAllDecayingParticles[loc_i].begin(), locAllDecayingParticles[loc_i].end());

		//Erase this vertex from future consideration
		locAllFullConstrainParticles.erase(locAllFullConstrainParticles.begin() + loc_i);
		locAllDecayingParticles.erase(locAllDecayingParticles.begin() + loc_i);
		locAllOnlyConstrainTimeParticles.erase(locAllOnlyConstrainTimeParticles.begin() + loc_i);
		locAllNoConstrainParticles.erase(locAllNoConstrainParticles.begin() + loc_i);
	}

	//all remaining vertex constraints now need decaying particles to help constrain.
	if(locAllFullConstrainParticles.empty() || !dLinkVerticesFlag)
	{
		if(dDebugLevel > 10)
			cout << "DKinFitUtils: Vertex constraints created." << endl;
		return locSortedConstraints; //either no remaining vertices, or cannot constrain them since linking via decaying particles is disabled
	}

	//loop over vertex-constraints-to-sort:
		//find which constraints decaying particles should be defined-by/constrained-to
		//find order in which constraints need to be constrained
	size_t locConstraintIndex = 0;
	bool locProgessMadeFlag = false;
	while(!locAllFullConstrainParticles.empty())
	{
		if(locConstraintIndex == locAllFullConstrainParticles.size())
		{
			if(!locProgessMadeFlag)
				break; //no progress made: cannot constrain remaining vertices
			locConstraintIndex = 0;
			locProgessMadeFlag = false;
			continue;
		}

		//find which decaying particles at this vertex have been previously defined
		set<DKinFitParticle*> locVertexDecayingParticles_Defined;
		set_intersection(locAllDecayingParticles[locConstraintIndex].begin(), locAllDecayingParticles[locConstraintIndex].end(),
                          locDefinedDecayingParticles.begin(), locDefinedDecayingParticles.end(),
                          inserter(locVertexDecayingParticles_Defined, locVertexDecayingParticles_Defined.begin()));

		//see if enough defined particles to constrain vertex
		if(locVertexDecayingParticles_Defined.size() + locAllFullConstrainParticles[locConstraintIndex].size() < 2)
		{
			++locConstraintIndex;
			continue; //nope
		}

		//find which decaying particles at this vertex have NOT been previously defined
		set<DKinFitParticle*> locVertexDecayingParticles_NotDefined;
		set_intersection(locAllDecayingParticles[locConstraintIndex].begin(), locAllDecayingParticles[locConstraintIndex].end(),
                          locDefinedDecayingParticles.begin(), locDefinedDecayingParticles.end(),
                          inserter(locVertexDecayingParticles_NotDefined, locVertexDecayingParticles_NotDefined.begin()));

		//Add decaying particles to appropriate sets
		locAllFullConstrainParticles[locConstraintIndex].insert(locVertexDecayingParticles_Defined.begin(), locVertexDecayingParticles_Defined.end());
		locAllNoConstrainParticles[locConstraintIndex].insert(locVertexDecayingParticles_NotDefined.begin(), locVertexDecayingParticles_NotDefined.end());

		//Create vertex/spacetime constraint, with decaying particles as no-constrain particles
		if(locSpacetimeFitFlag)
			locSortedConstraints.push_back(Make_SpacetimeConstraint(locAllFullConstrainParticles[locConstraintIndex], locAllOnlyConstrainTimeParticles[locConstraintIndex], locAllNoConstrainParticles[locConstraintIndex]));
		else //vertex only
		{
			locAllNoConstrainParticles[locConstraintIndex].insert(locAllOnlyConstrainTimeParticles[locConstraintIndex].begin(), locAllOnlyConstrainTimeParticles[locConstraintIndex].end());
			locSortedConstraints.push_back(Make_VertexConstraint(locAllFullConstrainParticles[locConstraintIndex], locAllNoConstrainParticles[locConstraintIndex]));
		}

		//The positions of these decaying particles are now defined: Can use to constrain vertices in later constraints
		locDefinedDecayingParticles.insert(locVertexDecayingParticles_NotDefined.begin(), locVertexDecayingParticles_NotDefined.end());

		//Erase this vertex from future consideration
		locAllFullConstrainParticles.erase(locAllFullConstrainParticles.begin() + locConstraintIndex);
		locAllDecayingParticles.erase(locAllDecayingParticles.begin() + locConstraintIndex);
		locAllOnlyConstrainTimeParticles.erase(locAllOnlyConstrainTimeParticles.begin() + locConstraintIndex);
		locAllNoConstrainParticles.erase(locAllNoConstrainParticles.begin() + locConstraintIndex);

		locProgessMadeFlag = true;
	}

	if(dDebugLevel > 10)
		cout << "DKinFitUtils: Vertex constraints created." << endl;

	return locSortedConstraints;
}

void DKinFitUtils::Group_VertexParticles(const set<DKinFitParticle*>& locVertexParticles, set<DKinFitParticle*>& locFullConstrainParticles, set<DKinFitParticle*>& locDecayingParticles, set<DKinFitParticle*>& locOnlyConstrainTimeParticles, set<DKinFitParticle*>& locNoConstrainParticles) const
{
	//Decaying: Only those that can conceivably be used to constrain: All unless dLinkVerticesFlag disabled (then none)
	locFullConstrainParticles.clear();
	locDecayingParticles.clear();
	locOnlyConstrainTimeParticles.clear();
	locNoConstrainParticles.clear();

	set<DKinFitParticle*>::const_iterator locIterator = locVertexParticles.begin();
	for(; locIterator != locVertexParticles.end(); ++locIterator)
	{
		DKinFitParticle* locParticle = *locIterator;
		DKinFitParticleType locKinFitParticleType = locParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
			locNoConstrainParticles.insert(locParticle);
		else if(locParticle->Get_KinFitParticleType() == d_BeamParticle)
		{
			if(Get_IncludeBeamlineInVertexFitFlag())
				locFullConstrainParticles.insert(locParticle);
			else
				locOnlyConstrainTimeParticles.insert(locParticle);
		}
		else if(locKinFitParticleType == d_DecayingParticle)
		{
			if(dLinkVerticesFlag)
				locDecayingParticles.insert(locParticle);
			else
				locNoConstrainParticles.insert(locParticle);
		}
		else if(locParticle->Get_Charge() == 0) //detected neutral
			locOnlyConstrainTimeParticles.insert(locParticle);
		else //detected charged
			locFullConstrainParticles.insert(locParticle);
	}
}

/**************************************************************** CREATE CONSTRAINTS ***************************************************************/

DKinFitConstraint_Mass* DKinFitUtils::Make_MassConstraint(DKinFitParticle* locDecayingParticle)
{
	if(locDecayingParticle->dKinFitParticleType != d_DecayingParticle)
	{
		cout << "ERROR: Wrong particle type in DKinFitUtils::Make_MassConstraint(). Returning NULL." << endl;
		return NULL;
	}

	//if constraint already exists for this particle, return it
	map<DKinFitParticle*, DKinFitConstraint_Mass*>::iterator locIterator = dMassConstraintMap.find(locDecayingParticle);
	if(locIterator != dMassConstraintMap.end())
		return locIterator->second;

	DKinFitConstraint_Mass* locConstraint = Get_KinFitConstraintMassResource();
	locConstraint->Set_DecayingParticle(locDecayingParticle);

	dMassConstraintMap[locDecayingParticle] = locConstraint;
	return locConstraint;
}

DKinFitConstraint_P4* DKinFitUtils::Make_P4Constraint(const set<DKinFitParticle*>& locInitialParticles, const set<DKinFitParticle*>& locFinalParticles)
{
	//if constraint already exists for this input, return it
	pair<set<DKinFitParticle*>, set<DKinFitParticle*> > locInputPair(locInitialParticles, locFinalParticles);
	map<pair<set<DKinFitParticle*>, set<DKinFitParticle*> >, DKinFitConstraint_P4*>::iterator locIterator = dP4ConstraintMap.find(locInputPair);
	if(locIterator != dP4ConstraintMap.end())
		return locIterator->second;

	DKinFitConstraint_P4* locConstraint = Get_KinFitConstraintP4Resource();
	locConstraint->Set_InitialParticles(locInitialParticles);
	locConstraint->Set_FinalParticles(locFinalParticles);

	dP4ConstraintMap[locInputPair] = locConstraint;
	return locConstraint;
}

DKinFitConstraint_Vertex* DKinFitUtils::Make_VertexConstraint(const set<DKinFitParticle*>& locFullConstrainParticles, const set<DKinFitParticle*>& locNoConstrainParticles, TVector3 locVertexGuess)
{
	//if constraint already exists for this input, return it
	pair<set<DKinFitParticle*>, set<DKinFitParticle*> > locInputPair(locFullConstrainParticles, locNoConstrainParticles);
	map<pair<set<DKinFitParticle*>, set<DKinFitParticle*> >, DKinFitConstraint_Vertex*>::iterator locIterator = dVertexConstraintMap.find(locInputPair);
	if(locIterator != dVertexConstraintMap.end())
		return locIterator->second;

	DKinFitConstraint_Vertex* locConstraint = Get_KinFitConstraintVertexResource();
	locConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locConstraint->Set_InitVertexGuess(locVertexGuess);

	dVertexConstraintMap[locInputPair] = locConstraint;
	return locConstraint;
}

DKinFitConstraint_Spacetime* DKinFitUtils::Make_SpacetimeConstraint(const set<DKinFitParticle*>& locFullConstrainParticles, const set<DKinFitParticle*>& locOnlyConstrainTimeParticles, const set<DKinFitParticle*>& locNoConstrainParticles, TLorentzVector locSpacetimeGuess)
{
	cout << "ERROR: SPACETIME CONSTRAINTS ARE NOT SUPPORTED YET. RETURNING NULL FROM DKinFitUtils::Make_SpacetimeConstraint()." << endl;
	return NULL;

	//if constraint already exists for this input, return it
	DSpacetimeParticles locSpacetimeParticles(locFullConstrainParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles);
	map<DSpacetimeParticles, DKinFitConstraint_Spacetime*>::iterator locIterator = dSpacetimeConstraintMap.find(locSpacetimeParticles);
	if(locIterator != dSpacetimeConstraintMap.end())
		return locIterator->second;

	DKinFitConstraint_Spacetime* locConstraint = Get_KinFitConstraintSpacetimeResource();
	locConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locConstraint->Set_OnlyConstrainTimeParticles(locOnlyConstrainTimeParticles);
	locConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locConstraint->Set_InitVertexGuess(locSpacetimeGuess.Vect());
	locConstraint->Set_InitTimeGuess(locSpacetimeGuess.T());

	dSpacetimeConstraintMap[locSpacetimeParticles] = locConstraint;
	return locConstraint;
}

/********************************************************** CLONE PARTICLES AND CONSTRAINTS ********************************************************/

TMatrixDSym* DKinFitUtils::Clone_MatrixDSym(const TMatrixDSym* locMatrix)
{
	if(locMatrix == NULL)
		return NULL;
	TMatrixDSym* locNewMatrix = Get_MatrixDSymResource();
	int locMatrixSize = locMatrix->GetNcols();
	locNewMatrix->ResizeTo(locMatrixSize, locMatrixSize);
	*locNewMatrix = *locMatrix;
	return locNewMatrix;
}

DKinFitParticle* DKinFitUtils::Clone_KinFitParticle(DKinFitParticle* locKinFitParticle)
{
	DKinFitParticle* locClonedKinFitParticle = Get_KinFitParticleResource();
	*locClonedKinFitParticle = *locKinFitParticle;
	dParticleMap_OutputToInput[locClonedKinFitParticle] = locKinFitParticle;

	if(dDebugLevel > 20)
		cout << "Cloned Particle: PID, input, output = " << locKinFitParticle->Get_PID() << ", " << locKinFitParticle << ", " << locClonedKinFitParticle << endl;

	//clone covariance matrix
	const TMatrixDSym* locCovarianceMatrix = locClonedKinFitParticle->Get_CovarianceMatrix();
	if(locCovarianceMatrix != NULL)
		locClonedKinFitParticle->Set_CovarianceMatrix(Clone_MatrixDSym(locCovarianceMatrix));

	return locClonedKinFitParticle;
}

DKinFitConstraint_P4* DKinFitUtils::Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_P4* locClonedConstraint = Get_KinFitConstraintP4Resource();
	*locClonedConstraint = *locConstraint;
//	dKinFitConstraintCloningOIMap[locClonedConstraint] = locConstraint;
	return locClonedConstraint;
}

DKinFitConstraint_Mass* DKinFitUtils::Clone_KinFitConstraint_Mass(const DKinFitConstraint_Mass* locConstraint)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_Mass* locClonedConstraint = Get_KinFitConstraintMassResource();
	*locClonedConstraint = *locConstraint;
//	dKinFitConstraintCloningOIMap[locClonedConstraint] = locConstraint;
	return locClonedConstraint;
}

DKinFitConstraint_Vertex* DKinFitUtils::Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_Vertex* locClonedConstraint = Get_KinFitConstraintVertexResource();
	*locClonedConstraint = *locConstraint;
//	dKinFitConstraintCloningOIMap[locClonedConstraint] = locConstraint;
	return locClonedConstraint;
}

DKinFitConstraint_Spacetime* DKinFitUtils::Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint)
{
	//to be called PRIOR to a fit
	DKinFitConstraint_Spacetime* locClonedConstraint = Get_KinFitConstraintSpacetimeResource();
	*locClonedConstraint = *locConstraint;
//	dKinFitConstraintCloningOIMap[locClonedConstraint] = locConstraint;
	return locClonedConstraint;
}

set<DKinFitParticle*> DKinFitUtils::Build_CloneParticleSet(const set<DKinFitParticle*>& locInputParticles, const map<DKinFitParticle*, DKinFitParticle*>& locCloneIOMap) const
{
	set<DKinFitParticle*> locCloneParticles;
	set<DKinFitParticle*>::iterator locParticleIterator = locInputParticles.begin();
	for(; locParticleIterator != locInputParticles.end(); ++locParticleIterator)
		locCloneParticles.insert(locCloneIOMap.find(*locParticleIterator)->second);
	return locCloneParticles;
}

set<DKinFitConstraint*> DKinFitUtils::Clone_ParticlesAndConstraints(const set<DKinFitConstraint*>& locInputConstraints)
{
	set<DKinFitConstraint*> locClonedConstraints;

	//Get all of the particles from the constraints (some particles may be listed in multiple constraints!)
		//This is why you can't clone the constraint particles one constraint at a time
	set<DKinFitParticle*> locAllParticles;
	set<DKinFitConstraint*>::const_iterator locConstraintIterator = locInputConstraints.begin();
	for(; locConstraintIterator != locInputConstraints.end(); ++locConstraintIterator)
	{
		set<DKinFitParticle*> locConstraintKinFitParticles = (*locConstraintIterator)->Get_AllParticles();
		locAllParticles.insert(locConstraintKinFitParticles.begin(), locConstraintKinFitParticles.end());

		//now, for those particles that may not directly be used in a constraint, but ARE used to define a decaying particle
		set<DKinFitParticle*>::iterator locParticleIterator = locConstraintKinFitParticles.begin();
		for(; locParticleIterator != locConstraintKinFitParticles.end(); ++locParticleIterator)
		{
			set<DKinFitParticle*> locFromAllButDecaying = (*locParticleIterator)->Get_FromAllButDecaying();
			locAllParticles.insert(locFromAllButDecaying.begin(), locFromAllButDecaying.end());
		}
	}

	//Clone all of the particles //keep track of clone IO for this fit
		//can't do as an overall class member, because one input may have several cloned outputs (multiple fits).
		//but for this fit, can track locally (here)
	map<DKinFitParticle*, DKinFitParticle*> locCloneIOMap; //for this fit
	set<DKinFitParticle*>::iterator locParticleIterator = locAllParticles.begin();
	for(; locParticleIterator != locAllParticles.end(); ++locParticleIterator)
		locCloneIOMap[*locParticleIterator] = Clone_KinFitParticle(*locParticleIterator);

	//Now, for all of the decaying cloned particles, go through and set new pointers for the from-initial and from-final state particles
	map<DKinFitParticle*, DKinFitParticle*>::iterator locCloneIterator = locCloneIOMap.begin();
	for(; locCloneIterator != locCloneIOMap.end(); ++locCloneIterator)
	{
		DKinFitParticle* locOutputParticle = locCloneIterator->second;
		if(locOutputParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //none

		//initial state
		set<DKinFitParticle*> locNewFromInitialState;
		set<DKinFitParticle*> locFromInitialState = locOutputParticle->Get_FromInitialState();
		set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
		for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
			locNewFromInitialState.insert(locCloneIOMap[*locParticleIterator]);
		locOutputParticle->Set_FromInitialState(locNewFromInitialState);

		//final state
		set<DKinFitParticle*> locNewFromFinalState;
		set<DKinFitParticle*> locFromFinalState = locOutputParticle->Get_FromFinalState();
		for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
			locNewFromFinalState.insert(locCloneIOMap[*locParticleIterator]);
		locOutputParticle->Set_FromFinalState(locNewFromFinalState);

	}

	//Clone the constraints, and then set the particles to the cloned particles
	for(locConstraintIterator = locInputConstraints.begin(); locConstraintIterator != locInputConstraints.end(); ++locConstraintIterator)
	{
		DKinFitConstraint_P4* locP4Constraint = dynamic_cast<DKinFitConstraint_P4*>(*locConstraintIterator);
		if(locP4Constraint != NULL)
		{
			DKinFitConstraint_P4* locClonedConstraint = Clone_KinFitConstraint_P4(locP4Constraint);
			locClonedConstraint->Set_InitialParticles(Build_CloneParticleSet(locClonedConstraint->Get_InitialParticles(), locCloneIOMap));
			locClonedConstraint->Set_FinalParticles(Build_CloneParticleSet(locClonedConstraint->Get_FinalParticles(), locCloneIOMap));
			locClonedConstraints.insert(locClonedConstraint);
			continue;
		}

		DKinFitConstraint_Mass* locMassConstraint = dynamic_cast<DKinFitConstraint_Mass*>(*locConstraintIterator);
		if(locMassConstraint != NULL)
		{
			DKinFitConstraint_Mass* locClonedConstraint = Clone_KinFitConstraint_Mass(locMassConstraint);
			locClonedConstraint->Set_DecayingParticle(locCloneIOMap.find(locClonedConstraint->Get_DecayingParticle())->second);
			locClonedConstraints.insert(locClonedConstraint);
			continue;
		}

		DKinFitConstraint_Vertex* locVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*locConstraintIterator);
		DKinFitConstraint_Spacetime* locSpacetimeConstraint = dynamic_cast<DKinFitConstraint_Spacetime*>(*locConstraintIterator);
		if((locVertexConstraint != NULL) && (locSpacetimeConstraint == NULL))
		{
			DKinFitConstraint_Vertex* locClonedConstraint = Clone_KinFitConstraint_Vertex(locVertexConstraint);
			locClonedConstraint->Set_FullConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_FullConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_NoConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_NoConstrainParticles(), locCloneIOMap));
			locClonedConstraints.insert(locClonedConstraint);
			continue;
		}

		if(locSpacetimeConstraint != NULL)
		{
			DKinFitConstraint_Spacetime* locClonedConstraint = Clone_KinFitConstraint_Spacetime(locSpacetimeConstraint);
			locClonedConstraint->Set_FullConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_FullConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_NoConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_NoConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_OnlyConstrainTimeParticles(Build_CloneParticleSet(locClonedConstraint->Get_OnlyConstrainTimeParticles(), locCloneIOMap));
			locClonedConstraints.insert(locClonedConstraint);
			continue;
		}
	}

	return locClonedConstraints;
}

/*************************************************************** VALIDATE CONSTRAINTS **************************************************************/

bool DKinFitUtils::Validate_Constraints(const set<DKinFitConstraint*>& locKinFitConstraints) const
{
	//Do independent of the kinematic fitter!
	//Empty for now. User can override in derived class
	return true;
}

/*************************************************************** CALCULATION ROUTINES **************************************************************/

TLorentzVector DKinFitUtils::Calc_DecayingP4_ByPosition(const DKinFitParticle* locKinFitParticle, bool locAtPositionFlag, bool locDontPropagateAtAllFlag) const
{
	//if input flag is true: return the value of the p4 at spot defined by locKinFitParticle->Get_Position() //else at the common vertex
		//useful for setting the momentum: locKinFitParticle->Set_Momentum()
	if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		return TLorentzVector();

	bool locP3DerivedAtProductionVertexFlag = !locKinFitParticle->Get_FromInitialState().empty(); //else decay vertex
	bool locP3DerivedAtPositionFlag = (locP3DerivedAtProductionVertexFlag == locKinFitParticle->Get_VertexP4AtProductionVertex());
	bool locDontPropagateDecayingP3Flag = (locP3DerivedAtPositionFlag == locAtPositionFlag);
	return Calc_DecayingP4(locKinFitParticle, locDontPropagateDecayingP3Flag, 1.0, locDontPropagateAtAllFlag);
}

TLorentzVector DKinFitUtils::Calc_DecayingP4_ByP3Derived(const DKinFitParticle* locKinFitParticle, bool locAtP3DerivedFlag, bool locDontPropagateAtAllFlag) const
{
	//if input flag is true: return the value of the p4 at the vertex where the p3-deriving particles are at
		//useful for doing mass constraints
	if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		return TLorentzVector();

	bool locDontPropagateDecayingP3Flag = locAtP3DerivedFlag;
	return Calc_DecayingP4(locKinFitParticle, locDontPropagateDecayingP3Flag, 1.0, locDontPropagateAtAllFlag);
}

TLorentzVector DKinFitUtils::Calc_DecayingP4_ByVertex(const DKinFitParticle* locKinFitParticle, bool locAtProductionVertexFlag, bool locDontPropagateAtAllFlag) const
{
	//if input flag is true: return the value of the p4 at the production vertex
		//else return it at the decay vertex
	if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		return TLorentzVector();

	bool locP3DerivedAtProductionVertexFlag = !locKinFitParticle->Get_FromInitialState().empty();
	bool locDontPropagateDecayingP3Flag = (locP3DerivedAtProductionVertexFlag == locAtProductionVertexFlag);
	return Calc_DecayingP4(locKinFitParticle, locDontPropagateDecayingP3Flag, 1.0, locDontPropagateAtAllFlag);
}

//Don't call this directly! Well you can, but it's a little confusing. Better to call the wrapper functions. 
TLorentzVector DKinFitUtils::Calc_DecayingP4(const DKinFitParticle* locKinFitParticle, bool locDontPropagateDecayingP3Flag, double locStateSignMultiplier, bool locDontPropagateAtAllFlag) const
{
	//locDontPropagateDecayingP3Flag: if true: don't propagate first decaying particle p3 from defined vertex to the other vertex
	//E, px, py, pz
	int locCharge = locKinFitParticle->Get_Charge();
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locBField = Get_IsBFieldNearBeamline() ? Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();

	//This section is calculated assuming that the p4 is NEEDED at the COMMON vertex
		//if not, need a factor of -1 on delta-x, and on the derivatives wrst the vertices
		//e.g. for detected particles, needed p4 is at the common vertex, but it is defined at some other point on its trajectory
	//HOWEVER, for decaying particles, this MAY NOT be true: it MAY be defined at the position where it is needed
		//Decaying particles technically have two common vertices: where it is produced and where it decays
		//So here, the DEFINED vertex is where its position is defined by the other tracks, 
		//and its COMMON vertex is the vertex where it is used to constrain another vertex
	//e.g. g, p -> K+, K+, Xi-    Xi- -> pi-, Lambda    Lambda -> p, pi-
		//Assuming standard constraint setup: p4 constraint is initial step, mass constraints by invariant mass, Xi- vertex defined by K's
		//For Xi-, the p3 is defined at its decay vertex (by decay products)
		//And the Xi- DEFINED vertex is at its production vertex (from kaons)
		//But the p3 is NEEDED at the production vertex, which is where it's DEFINED
		//Thus we need a factor of -1
	bool locNeedP4AtProductionVertex = locKinFitParticle->Get_FromInitialState().empty(); //true if defined by decay products; else by missing mass
	double locVertexSignMultiplier = (locNeedP4AtProductionVertex == locKinFitParticle->Get_VertexP4AtProductionVertex()) ? -1.0 : 1.0;
	TVector3 locDeltaX = locVertexSignMultiplier*(locCommonVertex - locPosition); //vector points in the OPPOSITE direction of the momentum

	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	bool locCommonVertexFitFlag = locKinFitParticle->Get_FitCommonVertexFlag();
	bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();

	TLorentzVector locP4Sum;
	if(locKinFitParticleType != d_DecayingParticle)
		locP4Sum += locStateSignMultiplier*locP4; //all but enclosed decaying particle: will instead get p4 from decay products (may still need to propagate it below)

	if(dDebugLevel > 30)
		cout << "PID, sign, pxyzE = " << locKinFitParticle->Get_PID() << ", " << locStateSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(!locDontPropagateAtAllFlag && locCommonVertexFitFlag && locChargedBFieldFlag && ((locKinFitParticleType != d_DecayingParticle) || !locDontPropagateDecayingP3Flag))
	{
		//fitting vertex of charged track in magnetic field: momentum changes as function of vertex
		//decaying particles: p4 not directly used, deriving particles are: so must propagate if charged
		//if initial particle is a "detected" particle (actually a decaying particle treated as detected): still propagate vertex (assume p3/v3 defined at production vertex)

		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
		if(dDebugLevel > 30)
			cout << "propagate pxyz by: " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.X() << ", " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.Y() << ", " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.Z() << endl;

		locP4Sum.SetVect(locP4Sum.Vect() - locStateSignMultiplier*locA*locDeltaXCrossH);
	}

	if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_DecayingP4() Decaying Particle; PID = " << locKinFitParticle->Get_PID() << endl;

		//replace the decaying particle with the particles it's momentum is derived from
		//initial state
		set<DKinFitParticle*> locFromInitialState = locKinFitParticle->Get_FromInitialState();
		set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
		for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			locP4Sum += Calc_DecayingP4(*locParticleIterator, false, locStateSignMultiplier, locDontPropagateAtAllFlag); //decaying particle multiplier * 1.0
		}

		//final state
		set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			//If defined by invariant mass: add p4s of final state particles
			//If defined by missing mass: add p4s of init state, subtract final state
			locP4Sum += Calc_DecayingP4(*locParticleIterator, false, locNextStateSignMultiplier, locDontPropagateAtAllFlag);
		}
	}

	return locP4Sum;
}

bool DKinFitUtils::Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixDSym& locCovarianceMatrix) const
{
	// propagates the track info to the fit common vertex
		//returns false if nothing changed (info not propagated: e.g. missing particle), else returns true
	//assumes: that between the two points on the track (input measured point & kinfit common point): the b-field is constant and there is no eloss or multiple scattering 
	// this function acts in the following way on each particle, IN THIS ORDER OF EXECUTION:
		// if a common vertex was not fit, then the current results are returned
		// for the remaining types (including decaying particles):
			// the vertex is set to the common vertex
			// if the time was fit, then the time is set to the common time, else it is propagated to the common vertex
			// the path length is propagated
			// momentum:
				// if this is a neutral shower: momentum is redefined by the new vertex //already done by Update_ParticleParams()
				// if this is a neutral particle (either detected, beam, or decaying) or a charged particle without a b-field: the momentum is unchanged
				// if this is a charged particle in a b-field (either beam, detected, or decaying): the momentum is propagated to the common vertex
		//propagating the covariance matrix:
			//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
			//add common time to matrix: 11x11 or 9x9 (neutral shower): if kinfit just add in (no correlations to meas, corr to common v3), else transform
			//transform to 7x7: common v3 & common t are just copied to the measured spots
	//the output covariance matrix is 7x7, even if the particle represents a neutral shower (5x5)

	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if(!locKinFitParticle->Get_FitCommonVertexFlag())
		return false; // no distance over which to propagate

	if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
		return false; // particle properties already defined at the fit vertex

	locCovarianceMatrix = *(locKinFitParticle->Get_CovarianceMatrix());
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locBField = Get_IsBFieldNearBeamline() ? Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

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
		//note that if the common time is not kinfit, then the true uncertainty is overestimated: 
			//cannot be obtained without a kinematic fit (3 equations (xyz), one unknown (time))
	double locCommonTime;
	locCommonTParamIndex_TempMatrix = locCovarianceMatrix.GetNcols();
	if(locKinFitParticle->Get_FitCommonTimeFlag()) //spacetime was fit
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

bool DKinFitUtils::Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixDSym& locCovarianceMatrix, pair<double, double>& locPathLengthPair) const
{
	//locPathLengthPair: value, uncertainty
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if(!locKinFitParticle->Get_FitCommonVertexFlag())
		return false; // no distance over which to propagate

	if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
		return false; // particle properties already defined at the fit vertex

	int locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
	int locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locBField = Get_IsBFieldNearBeamline() ? Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
	TVector3 locH = locBField.Unit();

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
	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);
		double locPMag = locMomentum.Mag();
		locPathLengthPair.first = locDeltaXDotH*locPMag/locPDotH; //cos(theta_helix) = p.dot(h)/|p| = x.dot(h)/l (l = path length)

		TMatrixD locTransformationMatrix_PathLength(1, locTempMatrix.GetNcols());

//RESUME: DO I NEED TO BE CAREFUL ABOUT P3?? Probably
		TVector3 locDPathDP3 = (locDeltaXDotH/locPDotH)*((1.0/locPMag)*locMomentum - (locPMag/locPDotH)*locH);
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
		locPathLengthPair.second = sqrt(locTempMatrix(0, 0));
	}
	else // non-accelerating
	{
		locPathLengthPair.first = locDeltaX.Mag();

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
		locPathLengthPair.second = sqrt(locTempMatrix(0, 0));
	}

	return true;
}

void DKinFitUtils::Calc_DecayingParticleJacobian(const DKinFitParticle* locKinFitParticle, bool locDontPropagateDecayingP3Flag, double locStateSignMultiplier, int locNumEta, TMatrixD& locJacobian) const
{
	//locJacobian: matrix used to convert dV to the decaying particle covariance matrix: indices are px, py, pz, x, y, z, t
		//dimensions are: 7, (dNumXi + locNumEta);
	//uses defining-particles to calculate decaying particle information

	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	int locCharge = locKinFitParticle->Get_Charge();
	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locBField = Get_IsBFieldNearBeamline() ? Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	if(dDebugLevel > 50)
		cout << "jacobian: decay product: PID = " << locKinFitParticle->Get_PID() << endl;

	//This section is calculated assuming that the p4 is NEEDED at the COMMON vertex
		//if not, need a factor of -1 on delta-x, and on the derivatives wrst the vertices
		//e.g. for detected particles, needed p4 is at the common vertex, but it is defined at some other point on its trajectory
	//HOWEVER, for decaying particles, this MAY NOT be true: it MAY be defined at the position where it is needed
		//Decaying particles technically have two common vertices: where it is produced and where it decays
		//So here, the DEFINED vertex is where its position is defined by the other tracks, 
		//and its COMMON vertex is the vertex where it is used to constrain another vertex
	//e.g. g, p -> K+, K+, Xi-    Xi- -> pi-, Lambda    Lambda -> p, pi-
		//Assuming standard constraint setup: p4 constraint is initial step, mass constraints by invariant mass, Xi- vertex defined by K's
		//For Xi-, the p3 is defined at its decay vertex (by decay products)
		//And the Xi- DEFINED vertex is at its production vertex (from kaons)
		//But the p3 is NEEDED at the production vertex, which is where it's DEFINED
		//Thus we need a factor of -1
	bool locNeedP4AtProductionVertex = locKinFitParticle->Get_FromInitialState().empty(); //true if defined by decay products; else by missing mass
	double locVertexSignMultiplier = (locNeedP4AtProductionVertex == locKinFitParticle->Get_VertexP4AtProductionVertex()) ? -1.0 : 1.0;
	TVector3 locDeltaX = locVertexSignMultiplier*(locCommonVertex - locPosition); //vector points in the OPPOSITE direction of the momentum

	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	bool locCommonVertexFitFlag = locKinFitParticle->Get_FitCommonVertexFlag();
	bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
		locPxParamIndex += locNumEta;
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
		locVxParamIndex += locNumEta;
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex() + locNumEta;

	if(locKinFitParticleType == d_TargetParticle)
		return;
	else if(locChargedBFieldFlag && locCommonVertexFitFlag && (locKinFitParticleType != d_DecayingParticle))
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
		if(locChargedBFieldFlag && locKinFitParticle->Get_FitCommonVertexFlag() && !locDontPropagateDecayingP3Flag)
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

		//replace the decaying particle with the particles it's momentum is derived from
		//initial state
		set<DKinFitParticle*> locFromInitialState = locKinFitParticle->Get_FromInitialState();
		set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
		for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			Calc_DecayingParticleJacobian(*locParticleIterator, false, locStateSignMultiplier, locNumEta, locJacobian); //decaying particle multiplier * 1.0
		}

		//final state
		set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			//If defined by invariant mass: add p4s of final state particles
			//If defined by missing mass: add p4s of init state, subtract final state
			Calc_DecayingParticleJacobian(*locParticleIterator, false, locNextStateSignMultiplier, locNumEta, locJacobian);
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

const DKinFitChain* DKinFitUtils::Build_OutputKinFitChain(const DKinFitChain* locInputKinFitChain, set<DKinFitParticle*>& locKinFitOutputParticles)
{
	//First, build map of input -> output
	map<DKinFitParticle*, DKinFitParticle*> locInputToOutputParticleMap;
	set<DKinFitParticle*>::iterator locParticleIterator = locKinFitOutputParticles.begin();
	for(; locParticleIterator != locKinFitOutputParticles.end(); ++locParticleIterator)
		locInputToOutputParticleMap[dParticleMap_OutputToInput[*locParticleIterator]] = *locParticleIterator;

	DKinFitChain* locOutputKinFitChain = Get_KinFitChainResource();
	locOutputKinFitChain->Set_DefinedParticleStepIndex(locInputKinFitChain->Get_DefinedParticleStepIndex());
	locOutputKinFitChain->Set_IsInclusiveChannelFlag(locInputKinFitChain->Get_IsInclusiveChannelFlag());

	//loop over steps
	for(size_t loc_i = 0; loc_i < locInputKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
	{
		const DKinFitChainStep* locInputKinFitChainStep = locInputKinFitChain->Get_KinFitChainStep(loc_i);
		DKinFitChainStep* locOutputKinFitChainStep = Get_KinFitChainStepResource();

		locOutputKinFitChainStep->Set_InitialParticleDecayFromStepIndex(locInputKinFitChainStep->Get_InitialParticleDecayFromStepIndex());
		locOutputKinFitChainStep->Set_ConstrainDecayingMassFlag(locInputKinFitChainStep->Get_ConstrainDecayingMassFlag());

		set<DKinFitParticle*> locInitialParticles = locInputKinFitChainStep->Get_InitialParticles();
		for(locParticleIterator = locInitialParticles.begin(); locParticleIterator != locInitialParticles.end(); ++locParticleIterator)
		{
			if((*locParticleIterator)->Get_KinFitParticleType() == d_DecayingParticle)
				locOutputKinFitChain->Set_DecayStepIndex(*locParticleIterator, loc_i);

			map<DKinFitParticle*, DKinFitParticle*>::iterator locMapIterator = locInputToOutputParticleMap.find(*locParticleIterator);
			DKinFitParticle* locKinFitParticle = (locMapIterator != locInputToOutputParticleMap.end()) ? locMapIterator->second : *locParticleIterator;
			locOutputKinFitChainStep->Add_InitialParticle(locKinFitParticle);
		}

		set<DKinFitParticle*> locFinalParticles = locInputKinFitChainStep->Get_FinalParticles();
		for(locParticleIterator = locFinalParticles.begin(); locParticleIterator != locFinalParticles.end(); ++locParticleIterator)
		{
			map<DKinFitParticle*, DKinFitParticle*>::iterator locMapIterator = locInputToOutputParticleMap.find(*locParticleIterator);
			DKinFitParticle* locKinFitParticle = (locMapIterator != locInputToOutputParticleMap.end()) ? locMapIterator->second : *locParticleIterator;
			locOutputKinFitChainStep->Add_FinalParticle(locInputToOutputParticleMap[*locParticleIterator]);
		}

		locOutputKinFitChain->Add_KinFitChainStep(locOutputKinFitChainStep);
	}

	return locOutputKinFitChain;
}

