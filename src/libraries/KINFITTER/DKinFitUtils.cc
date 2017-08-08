#include "DKinFitUtils.h"

/*************************************************************** RESOURCE MANAGEMENT ***************************************************************/

DKinFitUtils::DKinFitUtils(void)
{
	dKinFitter = nullptr; //Is set by DKinFitter constructor
	dLinkVerticesFlag = true;
	dDebugLevel = 0;
	dUpdateCovarianceMatricesFlag = true;

	dResourcePool_TMatrixFSym = std::make_shared<DResourcePool<TMatrixFSym>>();
	dResourcePool_TMatrixFSym->Set_ControlParams(100, 20, 1000, 50000, 0);

	dResourcePool_KinFitParticle = std::make_shared<DResourcePool<DKinFitParticle>>();
	dResourcePool_KinFitParticle->Set_ControlParams(50, 20, 500, 5000, 0);

	dResourcePool_KinFitChainStep = std::make_shared<DResourcePool<DKinFitChainStep>>();
	dResourcePool_KinFitChainStep->Set_ControlParams(5, 5, 50, 100, 0);

	dResourcePool_KinFitChain = std::make_shared<DResourcePool<DKinFitChain>>();
	dResourcePool_KinFitChain->Set_ControlParams(5, 5, 50, 100, 0);

	dResourcePool_MassConstraint = std::make_shared<DResourcePool<DKinFitConstraint_Mass>>();
	dResourcePool_MassConstraint->Set_ControlParams(20, 20, 100, 300, 0);

	dResourcePool_P4Constraint = std::make_shared<DResourcePool<DKinFitConstraint_P4>>();
	dResourcePool_P4Constraint->Set_ControlParams(20, 20, 100, 300, 0);

	dResourcePool_VertexConstraint = std::make_shared<DResourcePool<DKinFitConstraint_Vertex>>();
	dResourcePool_VertexConstraint->Set_ControlParams(20, 20, 100, 300, 0);

	dResourcePool_SpacetimeConstraint = std::make_shared<DResourcePool<DKinFitConstraint_Spacetime>>();
	dResourcePool_SpacetimeConstraint->Set_ControlParams(20, 20, 100, 300, 0);
}

void DKinFitUtils::Reset_NewEvent(void)
{
	dParticleMap_OutputToInput.clear();
	dMassConstraintMap.clear();
	dP4ConstraintMap.clear();
	dVertexConstraintMap.clear();
	dSpacetimeConstraintMap.clear();
	Reset_NewFit();
}

shared_ptr<TMatrixFSym> DKinFitUtils::Get_SymMatrixResource(unsigned int locNumMatrixRows)
{
	auto locSymMatrix = dResourcePool_TMatrixFSym->Get_SharedResource();
	locSymMatrix->ResizeTo(locNumMatrixRows, locNumMatrixRows);
	return locSymMatrix;
}

/***************************************************************** CREATE PARTICLES ****************************************************************/

shared_ptr<DKinFitParticle> DKinFitUtils::Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex);

	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(locCovarianceMatrix);

	locKinFitParticle->Set_KinFitParticleType(d_BeamParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Beam particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils::Make_TargetParticle(int locPID, int locCharge, double locMass)
{
	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
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

shared_ptr<DKinFitParticle> DKinFitUtils::Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, double locPathLength, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 7) || (locCovarianceMatrix->GetNcols() != 7))
		return NULL; //is not 7x7

	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(locCharge);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_Momentum(locMomentum);
	locKinFitParticle->Set_CovarianceMatrix(locCovarianceMatrix);
	locKinFitParticle->Set_PathLength(locPathLength);
	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Detected particle created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils::Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix)
{
	if((locCovarianceMatrix->GetNrows() != 5) || (locCovarianceMatrix->GetNcols() != 5))
		return NULL; //is not 5x5

	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
	locKinFitParticle->Set_PID(locPID);
	locKinFitParticle->Set_Charge(0);
	locKinFitParticle->Set_Mass(locMass);
	locKinFitParticle->Set_IsNeutralShowerFlag(true);
	locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
	locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex); //will be updated with vertex guess
	locKinFitParticle->Set_ShowerEnergy(locShowerEnergy);
	locKinFitParticle->Set_CovarianceMatrix(locCovarianceMatrix);

	locKinFitParticle->Set_KinFitParticleType(d_DetectedParticle);

	if(dDebugLevel > 5)
	{
		cout << "DKinFitUtils: Detected shower created. Printing:" << endl;
		locKinFitParticle->Print_ParticleParams();
	}

	return locKinFitParticle;
}

shared_ptr<DKinFitParticle> DKinFitUtils::Make_MissingParticle(int locPID, int locCharge, double locMass)
{
	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
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

shared_ptr<DKinFitParticle> DKinFitUtils::Make_DecayingParticle(int locPID, int locCharge, double locMass, const set<shared_ptr<DKinFitParticle>>& locFromInitialState, const set<shared_ptr<DKinFitParticle>>& locFromFinalState)
{
	auto locKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
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

/**************************************************************** CREATE CONSTRAINTS ***************************************************************/

shared_ptr<DKinFitConstraint_Mass> DKinFitUtils::Make_MassConstraint(const shared_ptr<DKinFitParticle>& locDecayingParticle)
{
	if(locDecayingParticle->dKinFitParticleType != d_DecayingParticle)
	{
		cout << "ERROR: Wrong particle type in DKinFitUtils::Make_MassConstraint(). Returning NULL." << endl;
		return nullptr;
	}

	//if constraint already exists for this particle, return it
	auto locIterator = dMassConstraintMap.find(locDecayingParticle);
	if(locIterator != dMassConstraintMap.end())
		return locIterator->second;

	auto locConstraint = dResourcePool_MassConstraint->Get_SharedResource();
	locConstraint->Set_DecayingParticle(locDecayingParticle);

	dMassConstraintMap[locDecayingParticle] = locConstraint;
	return locConstraint;
}

shared_ptr<DKinFitConstraint_P4> DKinFitUtils::Make_P4Constraint(const set<shared_ptr<DKinFitParticle>>& locInitialParticles, const set<shared_ptr<DKinFitParticle>>& locFinalParticles)
{
	//if constraint already exists for this input, return it
	auto locInputPair = std::make_pair(locInitialParticles, locFinalParticles);
	auto locIterator = dP4ConstraintMap.find(locInputPair);
	if(locIterator != dP4ConstraintMap.end())
		return locIterator->second;

	auto locConstraint = dResourcePool_P4Constraint->Get_SharedResource();
	locConstraint->Set_InitialParticles(locInitialParticles);
	locConstraint->Set_FinalParticles(locFinalParticles);

	dP4ConstraintMap[locInputPair] = locConstraint;
	return locConstraint;
}

shared_ptr<DKinFitConstraint_Vertex> DKinFitUtils::Make_VertexConstraint(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles, const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles, TVector3 locVertexGuess)
{
	//if constraint already exists for this input, return it
	auto locInputPair = std::make_pair(locFullConstrainParticles, locNoConstrainParticles);
	auto locIterator = dVertexConstraintMap.find(locInputPair);
	if(locIterator != dVertexConstraintMap.end())
		return locIterator->second;

	auto locConstraint = dResourcePool_VertexConstraint->Get_SharedResource();
	locConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locConstraint->Set_InitVertexGuess(locVertexGuess);

	dVertexConstraintMap[locInputPair] = locConstraint;
	return locConstraint;
}

shared_ptr<DKinFitConstraint_Spacetime> DKinFitUtils::Make_SpacetimeConstraint(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles, const set<shared_ptr<DKinFitParticle>>& locOnlyConstrainTimeParticles, const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles, TLorentzVector locSpacetimeGuess)
{
	cout << "ERROR: SPACETIME CONSTRAINTS ARE NOT SUPPORTED YET. RETURNING NULL FROM DKinFitUtils::Make_SpacetimeConstraint()." << endl;
	return NULL;

	//if constraint already exists for this input, return it
	DSpacetimeParticles locSpacetimeParticles(locFullConstrainParticles, locOnlyConstrainTimeParticles, locNoConstrainParticles);
	auto locIterator = dSpacetimeConstraintMap.find(locSpacetimeParticles);
	if(locIterator != dSpacetimeConstraintMap.end())
		return locIterator->second;

	auto locConstraint = dResourcePool_SpacetimeConstraint->Get_SharedResource();
	locConstraint->Set_FullConstrainParticles(locFullConstrainParticles);
	locConstraint->Set_OnlyConstrainTimeParticles(locOnlyConstrainTimeParticles);
	locConstraint->Set_NoConstrainParticles(locNoConstrainParticles);
	locConstraint->Set_InitVertexGuess(locSpacetimeGuess.Vect());
	locConstraint->Set_InitTimeGuess(locSpacetimeGuess.T());

	dSpacetimeConstraintMap[locSpacetimeParticles] = locConstraint;
	return locConstraint;
}

/********************************************************** CLONE PARTICLES AND CONSTRAINTS ********************************************************/

shared_ptr<TMatrixFSym> DKinFitUtils::Clone_SymMatrix(const TMatrixFSym* locMatrix)
{
	if(locMatrix == NULL)
		return NULL;
	int locMatrixSize = locMatrix->GetNcols();
	auto locNewMatrix = Get_SymMatrixResource(locMatrixSize);
	*locNewMatrix = *locMatrix;
	return locNewMatrix;
}

shared_ptr<DKinFitParticle> DKinFitUtils::Clone_KinFitParticle(const shared_ptr<DKinFitParticle>& locKinFitParticle)
{
	auto locClonedKinFitParticle = dResourcePool_KinFitParticle->Get_SharedResource();
	*locClonedKinFitParticle = *locKinFitParticle;
	dParticleMap_OutputToInput[locClonedKinFitParticle] = locKinFitParticle;

	if(dDebugLevel > 20)
		cout << "Cloned Particle: PID, input, output = " << locKinFitParticle->Get_PID() << ", " << locKinFitParticle << ", " << locClonedKinFitParticle << endl;

	//clone covariance matrix
	auto locCovarianceMatrix = locClonedKinFitParticle->Get_CovarianceMatrix();
	if((locCovarianceMatrix != NULL) && dUpdateCovarianceMatricesFlag)
		locClonedKinFitParticle->Set_CovarianceMatrix(Clone_SymMatrix(locCovarianceMatrix.get()));

	return locClonedKinFitParticle;
}

shared_ptr<DKinFitConstraint_P4> DKinFitUtils::Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint)
{
	//to be called PRIOR to a fit
	auto locClonedConstraint = dResourcePool_P4Constraint->Get_SharedResource();
	*locClonedConstraint = *locConstraint;
	return locClonedConstraint;
}

shared_ptr<DKinFitConstraint_Mass> DKinFitUtils::Clone_KinFitConstraint_Mass(const DKinFitConstraint_Mass* locConstraint)
{
	//to be called PRIOR to a fit
	auto locClonedConstraint = dResourcePool_MassConstraint->Get_SharedResource();
	*locClonedConstraint = *locConstraint;
	return locClonedConstraint;
}

shared_ptr<DKinFitConstraint_Vertex> DKinFitUtils::Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint)
{
	//to be called PRIOR to a fit
	auto locClonedConstraint = dResourcePool_VertexConstraint->Get_SharedResource();
	*locClonedConstraint = *locConstraint;
	return locClonedConstraint;
}

shared_ptr<DKinFitConstraint_Spacetime> DKinFitUtils::Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint)
{
	//to be called PRIOR to a fit
	auto locClonedConstraint = dResourcePool_SpacetimeConstraint->Get_SharedResource();
	*locClonedConstraint = *locConstraint;
	return locClonedConstraint;
}

set<shared_ptr<DKinFitParticle>> DKinFitUtils::Build_CloneParticleSet(const set<shared_ptr<DKinFitParticle>>& locInputParticles, const map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitParticle>>& locCloneIOMap) const
{
	set<shared_ptr<DKinFitParticle>> locCloneParticles;
	for(auto& locParticle : locInputParticles)
		locCloneParticles.insert(locCloneIOMap.find(locParticle)->second);
	return locCloneParticles;
}

set<shared_ptr<DKinFitConstraint>> DKinFitUtils::Clone_ParticlesAndConstraints(const set<shared_ptr<DKinFitConstraint>>& locInputConstraints)
{
	set<shared_ptr<DKinFitConstraint>> locClonedConstraints;

	//Get all of the particles from the constraints (some particles may be listed in multiple constraints!)
		//This is why you can't clone the constraint particles one constraint at a time
	set<shared_ptr<DKinFitParticle>> locAllParticles;
	for(auto& locConstraint : locInputConstraints)
	{
		auto locConstraintKinFitParticles = locConstraint->Get_AllParticles();
		locAllParticles.insert(locConstraintKinFitParticles.begin(), locConstraintKinFitParticles.end());

		//now, for those particles that may not directly be used in a constraint, but ARE used to define a decaying particle
		for(auto& locParticle : locConstraintKinFitParticles)
		{
			auto locFromAllParticles = locParticle->Get_FromAllParticles();
			locAllParticles.insert(locFromAllParticles.begin(), locFromAllParticles.end());
		}
	}

	//Clone all of the particles //keep track of clone IO for this fit
		//can't do as an overall class member, because one input may have several cloned outputs (multiple fits).
		//but for this fit, can track locally (here)
	map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitParticle>> locCloneIOMap; //for this fit
	for(auto& locParticle : locAllParticles)
		locCloneIOMap[locParticle] = Clone_KinFitParticle(locParticle);

	//Now, for all of the decaying cloned particles, go through and set new pointers for the from-initial and from-final state particles
	for(auto& locClonePair : locCloneIOMap)
	{
		auto locOutputParticle = locClonePair.second;
		if(locOutputParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue; //none

		//initial state
		set<shared_ptr<DKinFitParticle>> locNewFromInitialState;
		auto locFromInitialState = locOutputParticle->Get_FromInitialState();
		for(auto& locParticle : locFromInitialState)
			locNewFromInitialState.insert(locCloneIOMap[locParticle]);
		locOutputParticle->Set_FromInitialState(locNewFromInitialState);

		//final state
		set<shared_ptr<DKinFitParticle>> locNewFromFinalState;
		auto locFromFinalState = locOutputParticle->Get_FromFinalState();
		for(auto& locParticle : locFromFinalState)
			locNewFromFinalState.insert(locCloneIOMap[locParticle]);
		locOutputParticle->Set_FromFinalState(locNewFromFinalState);
	}

	//Clone the constraints, and then set the particles to the cloned particles
	for(auto& locConstraint : locInputConstraints)
	{
		auto locP4Constraint = std::dynamic_pointer_cast<DKinFitConstraint_P4>(locConstraint);
		if(locP4Constraint != NULL)
		{
			auto locClonedConstraint = Clone_KinFitConstraint_P4(locP4Constraint.get());
			locClonedConstraint->Set_InitialParticles(Build_CloneParticleSet(locClonedConstraint->Get_InitialParticles(), locCloneIOMap));
			locClonedConstraint->Set_FinalParticles(Build_CloneParticleSet(locClonedConstraint->Get_FinalParticles(), locCloneIOMap));
			locClonedConstraints.insert(std::dynamic_pointer_cast<DKinFitConstraint>(locClonedConstraint));
			continue;
		}

		auto locMassConstraint = std::dynamic_pointer_cast<DKinFitConstraint_Mass>(locConstraint);
		if(locMassConstraint != NULL)
		{
			auto locClonedConstraint = Clone_KinFitConstraint_Mass(locMassConstraint.get());
			locClonedConstraint->Set_DecayingParticle(locCloneIOMap.find(locClonedConstraint->Get_DecayingParticle())->second);
			locClonedConstraints.insert(std::dynamic_pointer_cast<DKinFitConstraint>(locClonedConstraint));
			continue;
		}

		auto locVertexConstraint = std::dynamic_pointer_cast<DKinFitConstraint_Vertex>(locConstraint);
		auto locSpacetimeConstraint = std::dynamic_pointer_cast<DKinFitConstraint_Spacetime>(locConstraint);
		if((locVertexConstraint != NULL) && (locSpacetimeConstraint == NULL))
		{
			auto locClonedConstraint = Clone_KinFitConstraint_Vertex(locVertexConstraint.get());
			locClonedConstraint->Set_FullConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_FullConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_NoConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_NoConstrainParticles(), locCloneIOMap));
			locClonedConstraints.insert(std::dynamic_pointer_cast<DKinFitConstraint>(locClonedConstraint));
			continue;
		}

		if(locSpacetimeConstraint != NULL)
		{
			auto locClonedConstraint = Clone_KinFitConstraint_Spacetime(locSpacetimeConstraint.get());
			locClonedConstraint->Set_FullConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_FullConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_NoConstrainParticles(Build_CloneParticleSet(locClonedConstraint->Get_NoConstrainParticles(), locCloneIOMap));
			locClonedConstraint->Set_OnlyConstrainTimeParticles(Build_CloneParticleSet(locClonedConstraint->Get_OnlyConstrainTimeParticles(), locCloneIOMap));
			locClonedConstraints.insert(std::dynamic_pointer_cast<DKinFitConstraint>(locClonedConstraint));
			continue;
		}
	}

	return locClonedConstraints;
}

void DKinFitUtils::Recycle_LastFitMemory(set<shared_ptr<DKinFitConstraint>>& locKinFitConstraints)
{
	//Get all of the particles from the constraints (some particles may be listed in multiple constraints!)
	//This is why you can't clone the constraint particles one constraint at a time
	set<shared_ptr<DKinFitParticle>> locAllParticles;
	auto locConstraintIterator = locKinFitConstraints.begin();
	for(; locConstraintIterator != locKinFitConstraints.end(); ++locConstraintIterator)
	{
		auto locConstraintKinFitParticles = (*locConstraintIterator)->Get_AllParticles();
		locAllParticles.insert(locConstraintKinFitParticles.begin(), locConstraintKinFitParticles.end());

		//now, for those particles that may not directly be used in a constraint, but ARE used to define a decaying particle
		auto locParticleIterator = locConstraintKinFitParticles.begin();
		for(; locParticleIterator != locConstraintKinFitParticles.end(); ++locParticleIterator)
		{
			auto locFromAllParticles = (*locParticleIterator)->Get_FromAllParticles();
			locAllParticles.insert(locFromAllParticles.begin(), locFromAllParticles.end());
		}
	}

	//Remove them from the output-to-input clone map
	for(auto& locParticle : locAllParticles)
		dParticleMap_OutputToInput.erase(locParticle);

	locKinFitConstraints.clear();
}

/*************************************************************** VALIDATE CONSTRAINTS **************************************************************/

bool DKinFitUtils::Validate_Constraints(const set<shared_ptr<DKinFitConstraint>>& locKinFitConstraints) const
{
	//Do independent of the kinematic fitter!
	//Empty for now. User can override in derived class
	return true;
}

/*************************************************************** CALCULATION ROUTINES **************************************************************/

bool DKinFitUtils::Get_IsDecayingParticleDefinedByProducts(const DKinFitParticle* locKinFitParticle) const
{
	auto locFromInitState = locKinFitParticle->Get_FromInitialState();
	if(locFromInitState.empty())
		return true;
	if(locFromInitState.size() >= 2)
		return false;
	return ((*locFromInitState.begin())->Get_KinFitParticleType() != d_TargetParticle);
}

TLorentzVector DKinFitUtils::Calc_DecayingP4_ByPosition(const DKinFitParticle* locKinFitParticle, bool locAtPositionFlag, bool locDontPropagateAtAllFlag) const
{
	//if input flag is true: return the value of the p4 at spot defined by locKinFitParticle->Get_Position() //else at the common vertex
		//useful for setting the momentum: locKinFitParticle->Set_Momentum()
	if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
		return TLorentzVector();

	bool locP3DerivedAtProductionVertexFlag = !Get_IsDecayingParticleDefinedByProducts(locKinFitParticle); //else decay vertex
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

	bool locP3DerivedAtProductionVertexFlag = !Get_IsDecayingParticleDefinedByProducts(locKinFitParticle);
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
	bool locNeedP4AtProductionVertex = Get_IsDecayingParticleDefinedByProducts(locKinFitParticle); //true if defined by decay products; else by missing mass
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

	if(!locDontPropagateAtAllFlag && (locKinFitParticleType != d_MissingParticle) && (locKinFitParticleType != d_TargetParticle) && locCommonVertexFitFlag && locChargedBFieldFlag && ((locKinFitParticleType != d_DecayingParticle) || !locDontPropagateDecayingP3Flag))
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
		auto locFromInitialState = locKinFitParticle->Get_FromInitialState();
		for(auto& locParticle : locFromInitialState)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << locParticle->Get_PID() << endl;
			locP4Sum += Calc_DecayingP4(locParticle.get(), false, locStateSignMultiplier, locDontPropagateAtAllFlag); //decaying particle multiplier * 1.0
		}

		//final state
		auto locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(auto& locParticle : locFromFinalState)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << locParticle->Get_PID() << endl;
			//If defined by invariant mass: add p4s of final state particles
			//If defined by missing mass: add p4s of init state, subtract final state
			locP4Sum += Calc_DecayingP4(locParticle.get(), false, locNextStateSignMultiplier, locDontPropagateAtAllFlag);
		}
	}

	return locP4Sum;
}

bool DKinFitUtils::Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixFSym* locCovarianceMatrix) const
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

	auto locFitCovMatrix = locKinFitParticle->Get_CovarianceMatrix();
	if(dUpdateCovarianceMatricesFlag)
	{
		if(locFitCovMatrix != NULL)
			*locCovarianceMatrix = *locFitCovMatrix;
		else
			locCovarianceMatrix->Zero();
	}

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

	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();
	int locCommonTParamIndex = locKinFitParticle->Get_CommonTParamIndex();

	//FIRST, DO EVERYTHING BUT THE COVARIANCE MATRIX

	//common v3
	locSpacetimeVertex.SetVect(locCommonVertex);

	//common time
	double locCommonTime = 0.0;
	if(locKinFitParticle->Get_FitCommonTimeFlag()) //spacetime was fit
		locCommonTime = locKinFitParticle->Get_CommonTime();
	else if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locP4.Vect().Dot(locH);
		locCommonTime = locKinFitParticle->Get_Time() + locDeltaXDotH*locP4.E()/(29.9792458*locPDotH);
	}
	else if(!locNeutralShowerFlag) //non-accelerating, non-shower
	{
		double locDeltaXDotP = locDeltaX.Dot(locP4.Vect());
		locCommonTime = locKinFitParticle->Get_Time() + locDeltaXDotP*locP4.E()/(29.9792458*locP4.Vect().Mag2());
	}
	else //neutral shower
		locCommonTime = locKinFitParticle->Get_Time() - locDeltaX.Mag()*locP4.E()/(29.9792458*locP4.P());
	locSpacetimeVertex.SetT(locCommonTime);

	//p3
	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //charged & in b-field
		locMomentum = locP4.Vect() - locDeltaX.Cross(locA*locH);
	else //constant: either neutral or no b-field
		locMomentum = locP4.Vect();

	//if not updating the covariance matrix, skip to the path length and return
	if(!dUpdateCovarianceMatricesFlag)
		return Calc_PathLength(locKinFitParticle, locVXi, locCovarianceMatrix, locPathLengthPair);

	//UPDATE THE COVARIANCE MATRIX
	int locCommonVxParamIndex_TempMatrix, locCommonTParamIndex_TempMatrix;

	//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
	locCommonVxParamIndex_TempMatrix = locCovarianceMatrix->GetNcols();
	locCovarianceMatrix->ResizeTo(locCommonVxParamIndex_TempMatrix + 3, locCommonVxParamIndex_TempMatrix + 3);
	for(size_t loc_i = 0; loc_i < 3; ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < 3; ++loc_j)
			(*locCovarianceMatrix)(loc_i + locCommonVxParamIndex_TempMatrix, loc_j + locCommonVxParamIndex_TempMatrix) = (*locVXi)(locCommonVxParamIndex + loc_i, locCommonVxParamIndex + loc_j);
	}
	// done: no correlations between common vertex and measured params!!!

	//add common time to matrix: 11x11 or 9x9 (neutral shower): if kinfit just add in (no correlations to meas, corr to common v3), else transform
		//note that if the common time is not kinfit, then the true uncertainty is overestimated: 
			//cannot be obtained without a kinematic fit (3 equations (xyz), one unknown (time))
	locCommonTParamIndex_TempMatrix = locCovarianceMatrix->GetNcols();
	if(locKinFitParticle->Get_FitCommonTimeFlag()) //spacetime was fit
	{
		locCovarianceMatrix->ResizeTo(locCovarianceMatrix->GetNcols() + 1, locCovarianceMatrix->GetNcols() + 1);
		(*locCovarianceMatrix)(locCommonTParamIndex_TempMatrix, locCommonTParamIndex_TempMatrix) = (*locVXi)(locCommonTParamIndex, locCommonTParamIndex);
		for(size_t loc_i = 0; loc_i < 3; ++loc_i) //correlations to common v3
		{
			(*locCovarianceMatrix)(locCommonTParamIndex_TempMatrix, locCommonVxParamIndex_TempMatrix + loc_i) = (*locVXi)(locCommonTParamIndex, locCommonVxParamIndex + loc_i);
			(*locCovarianceMatrix)(locCommonVxParamIndex_TempMatrix + loc_i, locCommonTParamIndex_TempMatrix) = (*locVXi)(locCommonVxParamIndex + loc_i, locCommonTParamIndex);
		}
	}
	else if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locP4.Vect().Dot(locH);

		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix->GetNcols() + 1, locCovarianceMatrix->GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix->GetNcols(); ++loc_i)
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

		locCovarianceMatrix->Similarity(locTransformationMatrix_CommonTime);
	}
	else if(!locNeutralShowerFlag) //non-accelerating, non-shower
	{
		double locDeltaXDotP = locDeltaX.Dot(locP4.Vect());

		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix->GetNcols() + 1, locCovarianceMatrix->GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix->GetNcols(); ++loc_i)
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

		locCovarianceMatrix->Similarity(locTransformationMatrix_CommonTime);
	}
	else //neutral shower
	{
		TMatrixD locTransformationMatrix_CommonTime(locCovarianceMatrix->GetNcols() + 1, locCovarianceMatrix->GetNcols());
		for(unsigned int loc_i = 0; int(loc_i) < locCovarianceMatrix->GetNcols(); ++loc_i)
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

		locCovarianceMatrix->Similarity(locTransformationMatrix_CommonTime);
	}

	//transform to 7x7: common v3 & common t are just copied to the measured spots; p is propagated if in bfield, else is copied
	TMatrixD locTransformationMatrix_Propagation(7, locCovarianceMatrix->GetNcols());
	//p3
	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //charged & in b-field
	{
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
		for(unsigned int loc_i = 0; loc_i < 3; ++loc_i)
			locTransformationMatrix_Propagation(loc_i, locCovMatrixPxParamIndex + loc_i) = 1.0;
	}

	//v3
	for(unsigned int loc_i = 0; loc_i < 3; ++loc_i)
		locTransformationMatrix_Propagation(3 + loc_i, locCommonVxParamIndex_TempMatrix + loc_i) = 1.0;

	//t
	locTransformationMatrix_Propagation(6, locCommonTParamIndex_TempMatrix) = 1.0;

	//transform!!
	locCovarianceMatrix->Similarity(locTransformationMatrix_Propagation); //FINALLY!!!

	//now calculate the path length
	return Calc_PathLength(locKinFitParticle, locVXi, locCovarianceMatrix, locPathLengthPair);
}

bool DKinFitUtils::Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixFSym* locCovarianceMatrix, pair<double, double>& locPathLengthPair) const
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

	//First, compute the path length
	if((locCharge != 0) && Get_IsBFieldNearBeamline()) //in b-field & charged
	{
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);
		double locPMag = locMomentum.Mag();
		locPathLengthPair.first = locDeltaXDotH*locPMag/locPDotH; //cos(theta_helix) = p.dot(h)/|p| = x.dot(h)/l (l = path length)
	}
	else // non-accelerating
		locPathLengthPair.first = locDeltaX.Mag();

	//if not updating the errors, set the error to zero
	if((locCovarianceMatrix == nullptr) || !dUpdateCovarianceMatricesFlag)
	{
		locPathLengthPair.second = 0.0;
		return true;
	}

	//now compute the uncertainty
	//add common v3 to matrix: 10x10 or 8x8 (neutral shower)
	TMatrixFSym locTempMatrix(*locCovarianceMatrix);

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

		TMatrixD locTransformationMatrix_PathLength(1, locTempMatrix.GetNcols());

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
	}
	else // non-accelerating
	{
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
	}
	locPathLengthPair.second = sqrt(locTempMatrix(0, 0));

	return true;
}

void DKinFitUtils::Calc_DecayingParticleJacobian(const DKinFitParticle* locKinFitParticle, bool locDontPropagateDecayingP3Flag, double locStateSignMultiplier, int locNumEta, const map<const DKinFitParticle*, int>& locAdditionalPxParamIndices, TMatrixD& locJacobian) const
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
	bool locNeedP4AtProductionVertex = Get_IsDecayingParticleDefinedByProducts(locKinFitParticle); //true if defined by decay products; else by missing mass
	double locVertexSignMultiplier = (locNeedP4AtProductionVertex == locKinFitParticle->Get_VertexP4AtProductionVertex()) ? -1.0 : 1.0;
	TVector3 locDeltaX = locVertexSignMultiplier*(locCommonVertex - locPosition); //vector points in the OPPOSITE direction of the momentum

	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	bool locCommonVertexFitFlag = locKinFitParticle->Get_FitCommonVertexFlag();
	bool locChargedBFieldFlag = (locCharge != 0) && Get_IsBFieldNearBeamline();
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	if(((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle)) && (locPxParamIndex >= 0))
		locPxParamIndex += locNumEta;
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	if(((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle)) && (locVxParamIndex >= 0))
		locVxParamIndex += locNumEta;
	if(locPxParamIndex < 0)
	{
		//for particles not included in the fit matrices
		auto locPxParamIterator = locAdditionalPxParamIndices.find(locKinFitParticle);
		if(locPxParamIterator != locAdditionalPxParamIndices.end())
			locPxParamIndex = locPxParamIterator->second;
	}
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
		auto locFromInitialState = locKinFitParticle->Get_FromInitialState();
		for(auto& locParticle : locFromInitialState)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << locParticle->Get_PID() << endl;
			Calc_DecayingParticleJacobian(locParticle.get(), false, locStateSignMultiplier, locNumEta, locAdditionalPxParamIndices, locJacobian); //decaying particle multiplier * 1.0
		}

		//final state
		auto locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(auto& locParticle : locFromFinalState)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << locParticle->Get_PID() << endl;
			//If defined by invariant mass: add p4s of final state particles
			//If defined by missing mass: add p4s of init state, subtract final state
			Calc_DecayingParticleJacobian(locParticle.get(), false, locNextStateSignMultiplier, locNumEta, locAdditionalPxParamIndices, locJacobian);
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

shared_ptr<const DKinFitChain> DKinFitUtils::Build_OutputKinFitChain(const shared_ptr<const DKinFitChain>& locInputKinFitChain, set<shared_ptr<DKinFitParticle>>& locKinFitOutputParticles)
{
	if(dDebugLevel > 20)
	{
		cout << "DKinFitUtils::Build_OutputKinFitChain(): Printing input chain." << endl;
		locInputKinFitChain->Print_InfoToScreen();
	}

	//First, build map of input -> output
	map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitParticle>> locInputToOutputParticleMap;
	for(auto& locParticle : locKinFitOutputParticles)
		locInputToOutputParticleMap[dParticleMap_OutputToInput[locParticle]] = locParticle;

	auto locOutputKinFitChain = dResourcePool_KinFitChain->Get_SharedResource();
	locOutputKinFitChain->Set_DefinedParticleStepIndex(locInputKinFitChain->Get_DefinedParticleStepIndex());
	locOutputKinFitChain->Set_IsInclusiveChannelFlag(locInputKinFitChain->Get_IsInclusiveChannelFlag());

	//loop over steps
	for(size_t loc_i = 0; loc_i < locInputKinFitChain->Get_NumKinFitChainSteps(); ++loc_i)
	{
		auto locInputKinFitChainStep = locInputKinFitChain->Get_KinFitChainStep(loc_i);
		auto locOutputKinFitChainStep = dResourcePool_KinFitChainStep->Get_SharedResource();

		locOutputKinFitChainStep->Set_InitialParticleDecayFromStepIndex(locInputKinFitChainStep->Get_InitialParticleDecayFromStepIndex());
		locOutputKinFitChainStep->Set_ConstrainDecayingMassFlag(locInputKinFitChainStep->Get_ConstrainDecayingMassFlag());

		auto locInitialParticles = locInputKinFitChainStep->Get_InitialParticles();
		for(auto locKinFitParticle : locInitialParticles)
		{
			if(locKinFitParticle == nullptr)
				continue;
			auto locMapIterator = locInputToOutputParticleMap.find(locKinFitParticle);
			if(locMapIterator != locInputToOutputParticleMap.end())
				locKinFitParticle = locMapIterator->second;
			locOutputKinFitChainStep->Add_InitialParticle(locKinFitParticle);
			if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle)
				locOutputKinFitChain->Set_DecayStepIndex(locKinFitParticle, loc_i);
		}

		auto locFinalParticles = locInputKinFitChainStep->Get_FinalParticles();
		for(auto locKinFitParticle : locFinalParticles)
		{
			if(locKinFitParticle == nullptr)
				continue;
			auto locMapIterator = locInputToOutputParticleMap.find(locKinFitParticle);
			if(locMapIterator != locInputToOutputParticleMap.end())
				locKinFitParticle = locMapIterator->second;
			locOutputKinFitChainStep->Add_FinalParticle(locKinFitParticle);
		}

		locOutputKinFitChain->Add_KinFitChainStep(locOutputKinFitChainStep);
	}

	if(dDebugLevel > 20)
	{
		cout << "DKinFitUtils::Build_OutputKinFitChain(): Printing output chain." << endl;
		locOutputKinFitChain->Print_InfoToScreen();
	}

	return locOutputKinFitChain;
}
