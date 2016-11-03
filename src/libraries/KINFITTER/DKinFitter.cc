#include "DKinFitter.h"
#include "DKinFitUtils.h"

//MAGNETIC FIELD
	//If the particle is charged, the magnetic field is taken into account if a common vertex is simultaneously kinematically fit
		//(magnetic field assumed to be constant over the propagation range)
		//If a common vertex is not defined, the momentum is assumed to be constant
	//Energy loss between the common vertex and the particle position is ignored

//NEUTRAL SHOWERS
	//For input neutral showers, the momentum is defined-by and updated-with the common vertex (see "Common Vertex Notes")
		//if a common vertex is not simultaneously kinematically fit for a neutral shower, the fit will fail
		//To include a neutral shower in a p4 fit without a vertex fit, create a neutral particle from it 
			//with a full 7x7 covaraince matrix and input it instead of the neutral shower
	//Massive neutral showers (e.g. neutron) cannot be used in vertex constraints: only spacetime constraints. However, photons can. 
		//This is because their momentum is defined by the vertex time

//P4 CONSTRAINTS:
	//Don't do if studying an inclusive channel.
	//Pick the step of the decay chain that contains the missing or open-ended-decaying particle, to make sure it is constrained
		//If no missing or open-ended decaying particle, pick the initial step, to make sure the beam information is used

//MASS CONSTRAINTS:
	//Cannot use to constrain the mass of a missing particle: use a p4 constraint instead
	//Mass of decaying particle constrained by either missing mass or invariant mass
		//If no missing particles, can constrain either way
			//By default, prefer using the invariant mass. This is to ensure consistency when combining with a p4 constraint.
			//Failure to do so could lead to uninvertible matrices, and a failed kinematic fit. 
		//If this particle has a missing decay product (inclusive or otherwise): Constrain by missing mass
		//If computation by missing mass would have a missing particle (inclusive or otherwise): Constrain by invariant mass
	//However, the way it is calculated is not defined in these constraints; it is defined when creating the decaying DKinFitParticle
		//So be sure to create them correctly!

//VERTEX CONSTRAINTS:
	//You can include the beam particle in the common vertex fit, but only do so if it has a non-zero error matrix.
	//Neutral and missing particles included in the constraint will not be used to constrain the vertex, but will be set with the fit common vertex
		//This is necessary if you want the neutral particle momentum (from an input shower) to change with the reconstructed vertex
	//decaying particles should only be used to constrain a fit if the position is defined in another vertex constraint
	//Massive neutral showers (e.g. neutron) cannot be used in vertex constraints: only spacetime constraints. However, photons can. 
		//This is because their momentum is defined by the vertex time

//SPACETIME CONSTRAINTS:
	//THESE ARE CURRENTLY DISABLED
	//It is not possible to fit a common time without simultaneously fitting a common vertex.
		//Requiring a common time at a non-common vertex has no meaning, and the fitter is not setup to input a pre-fit common vertex with uncertainties.
	//You can include the RF time in the fit by calling the Make_BeamParticle function with the RF information
	//Missing particles included in the constraint will not be used to constrain the spacetime, but will be set with the fit common spacetime
	//decaying particles should only be used to constrain a fit if the spacetime is defined in another vertex constraint

//OBJECTS & MEMORY:
	//DKinFitParticle and DKinFitConstraint objects have private constructors: User cannot create them directly 
		//Instead, they are managed by DKinFitUtils
	//User creates DKinFitParticles, adds them to the constraints
	//When the fit starts, INTERNALLY ONLY, the particles are cloned and the constraints are cloned to contain the clones
	//User can resuse particles & constraints objects

/**************************************************************** CONSTRUCTOR AND RESET ****************************************************************/

DKinFitter::DKinFitter(DKinFitUtils* locKinFitUtils) : dKinFitUtils(locKinFitUtils)
{
	dDebugLevel = 0;

	dMaxNumIterations = 20;
	dConvergenceChiSqDiff = 0.001;
	dConvergenceChiSqDiff_LastResort = 0.005;

	dKinFitUtils->dKinFitter = this;
	Reset_NewEvent();
}

void DKinFitter::Reset_NewEvent(void)
{
	dKinFitUtils->Reset_NewEvent();
	Reset_NewFit();
}

void DKinFitter::Reset_NewFit(void)
{
	dKinFitUtils->Reset_NewFit();

	dKinFitStatus = d_KinFitSuccessful;

	dKinFitConstraints.clear();
	dKinFitParticles.clear();
	dParticleConstraintMap.clear();
	dParticleConstraintMap_Direct.clear();

	dNumXi = 0;
	dNumEta = 0;
	dNumF = 0;

	dChiSq = 0.0;
	dNDF = 0;
	dConfidenceLevel = 0.0;
	dPulls.clear();

	dV = dKinFitUtils->Get_LargeMatrixDSymResource();
	dVXi = dKinFitUtils->Get_LargeMatrixDSymResource();
	dVEta = dKinFitUtils->Get_LargeMatrixDSymResource();
}

/********************************************************************** UTILITIES **********************************************************************/

void DKinFitter::Set_DebugLevel(int locDebugLevel)
{
	dDebugLevel = locDebugLevel;
	dKinFitUtils->Set_DebugLevel(dDebugLevel);
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

bool DKinFitter::Get_IsVertexConstrained(DKinFitParticle* locKinFitParticle) const
{
	set<DKinFitConstraint_Vertex*> locConstraints = Get_Constraints<DKinFitConstraint_Vertex>(locKinFitParticle);
	set<DKinFitConstraint_Vertex*>::iterator locSetIterator = locConstraints.begin();
	for(; locSetIterator != locConstraints.end(); ++locSetIterator)
	{
		set<DKinFitParticle*> locConstrainedParticles = (*locSetIterator)->Get_AllConstrainedParticles();
		if(locConstrainedParticles.find(locKinFitParticle) != locConstrainedParticles.end())
			return true;
	}
	return false;
}

bool DKinFitter::Get_IsTimeConstrained(DKinFitParticle* locKinFitParticle) const
{
	set<DKinFitConstraint_Spacetime*> locConstraints = Get_Constraints<DKinFitConstraint_Spacetime>(locKinFitParticle);
	set<DKinFitConstraint_Spacetime*>::iterator locSetIterator = locConstraints.begin();
	for(; locSetIterator != locConstraints.end(); ++locSetIterator)
	{
		set<DKinFitParticle*> locConstrainedParticles = (*locSetIterator)->Get_AllConstrainedParticles();
		if(locConstrainedParticles.find(locKinFitParticle) != locConstrainedParticles.end())
			return true;
	}
	return false;
}

/******************************************************************* INITIALIZATION ********************************************************************/

void DKinFitter::Prepare_ConstraintsAndParticles(void)
{
	//Print constraint info
	if(dDebugLevel > 5)
	{
		cout << "DKinFitter: Pre-fit constraint info:" << endl;
		set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
		for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
			(*locConstraintIterator)->Print_ConstraintInfo();
	}

	//clone constraints & particles
	dKinFitConstraints = dKinFitUtils->Clone_ParticlesAndConstraints(dKinFitConstraints);

	//Print cloned constraint info
	if(dDebugLevel > 10)
	{
		cout << "DKinFitter: Cloned (Fit) constraint info:" << endl;
		set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
		for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
			(*locConstraintIterator)->Print_ConstraintInfo();
	}

	//Build dKinFitParticles, dParticleConstraintMap & _Direct
	set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
	for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
	{
		set<DKinFitParticle*> locKinFitParticles = (*locConstraintIterator)->Get_AllParticles();
		dKinFitParticles.insert(locKinFitParticles.begin(), locKinFitParticles.end());

		//loop over the particles
		set<DKinFitParticle*>::iterator locParticleIterator = locKinFitParticles.begin();
		for(; locParticleIterator != locKinFitParticles.end(); ++locParticleIterator)
		{
			dParticleConstraintMap[*locParticleIterator].insert(*locConstraintIterator);
			dParticleConstraintMap_Direct[*locParticleIterator].insert(*locConstraintIterator);

			if((*locParticleIterator)->Get_KinFitParticleType() != d_DecayingParticle)
				continue;

			//now, check for particles that may not directly be used in this constraint, but ARE used to define this decaying particle
			
			//this will not be the case if the particle is in a vertex constraint, but only as a no-constrain particle
			DKinFitConstraint_Vertex* locVertexConstraint = dynamic_cast<DKinFitConstraint_Vertex*>(*locConstraintIterator);
			if(locVertexConstraint != NULL)
			{
				set<DKinFitParticle*> locNoConstrainParticles = locVertexConstraint->Get_NoConstrainParticles();
				if(locNoConstrainParticles.find(*locParticleIterator) != locNoConstrainParticles.end())
					continue; //is not used to constrain, so neither is it's p4-defining particles (or at least, not here)
			}

			//this decaying particle is directly used in the constraint.  all p4-defining particles are then also used (but indirectly)
			set<DKinFitParticle*> locFromAllParticles = (*locParticleIterator)->Get_FromAllParticles();
			dKinFitParticles.insert(locFromAllParticles.begin(), locFromAllParticles.end());

			//dParticleConstraintMap
			set<DKinFitParticle*>::iterator locDerivingParticleIterator = locFromAllParticles.begin();
			for(; locDerivingParticleIterator != locFromAllParticles.end(); ++locDerivingParticleIterator)
				dParticleConstraintMap[*locDerivingParticleIterator].insert(*locConstraintIterator);
		}
	}

	//Set dVertexP4AtProductionVertex on decaying particles
	//Loop over vertex constraints, searching for decaying particles where their positions are defined
	set<DKinFitConstraint_Vertex*> locVertexConstraints = Get_Constraints<DKinFitConstraint_Vertex>();
	set<DKinFitConstraint_Vertex*>::iterator locVertexIterator = locVertexConstraints.begin();
	for(; locVertexIterator != locVertexConstraints.end(); ++locVertexIterator)
	{
		set<DKinFitParticle*> locAllConstraintParticles = (*locVertexIterator)->Get_AllParticles();
		set<DKinFitParticle*> locNoConstrainParticles = (*locVertexIterator)->Get_NoConstrainParticles();
		set<DKinFitParticle*>::iterator locParticleIterator = locNoConstrainParticles.begin();
		for(; locParticleIterator != locNoConstrainParticles.end(); ++locParticleIterator)
		{
			DKinFitParticle* locKinFitParticle = *locParticleIterator;
			if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
				continue;

			//We have the decaying particle where it's position is defined (no-constrain)
			//Now we must determine whether this is the production or decay vertex

			//See which from-final-state particles are found in the vertex constraint
			set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
			set<DKinFitParticle*> locFinalStateParticlesAtVertex;
			set_intersection(locAllConstraintParticles.begin(), locAllConstraintParticles.end(), locFromFinalState.begin(), 
				locFromFinalState.end(), inserter(locFinalStateParticlesAtVertex, locFinalStateParticlesAtVertex.begin()));

			//Set vertex flag based on whether any particles found or not
			bool locIsProductionVertex = false;
			if(locKinFitParticle->Get_FromInitialState().empty()) //p4 defined by invariant mass (decay products)
				locIsProductionVertex = locFinalStateParticlesAtVertex.empty() ? true : false; //if no overlap: not at decay vertex: production
			else //p4 defined by missing mass
				locIsProductionVertex = locFinalStateParticlesAtVertex.empty() ? false : true; //if no overlap: not at production vertex: decay
			locKinFitParticle->Set_VertexP4AtProductionVertex(locIsProductionVertex);
		}
	}

	//Initialize un-set particle data

	//Common vertices
	for(locVertexIterator = locVertexConstraints.begin(); locVertexIterator != locVertexConstraints.end(); ++locVertexIterator)
		(*locVertexIterator)->Set_CommonVertex((*locVertexIterator)->Get_InitVertexGuess());

	//Common times (init vertex guess already set above)
	set<DKinFitConstraint_Spacetime*> locSpacetimeConstraints = Get_Constraints<DKinFitConstraint_Spacetime>();
	set<DKinFitConstraint_Spacetime*>::iterator locSpacetimeIterator = locSpacetimeConstraints.begin();
	for(; locSpacetimeIterator != locSpacetimeConstraints.end(); ++locSpacetimeIterator)
		(*locSpacetimeIterator)->Set_CommonTime((*locSpacetimeIterator)->Get_InitTimeGuess());

	//neutral shower p3
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locParticle = *locParticleIterator;
		if(!locParticle->Get_IsNeutralShowerFlag())
			continue; //only do for neutral showers

		double locE = locParticle->Get_ShowerEnergy();
		double locMass = locParticle->Get_Mass();
		double locPMag = sqrt(locE*locE - locMass*locMass);
		TVector3 locMomentum = locParticle->Get_Position() - locParticle->Get_CommonVertex();
		locMomentum.SetMag(locPMag);
		locParticle->Set_Momentum(locMomentum);
	}

	//missing/open-ended-decaying p3
	set<DKinFitConstraint_P4*> locP4Constraints = Get_Constraints<DKinFitConstraint_P4>();
	if(!locP4Constraints.empty())
	{
		DKinFitConstraint_P4* locP4Constraint = *(locP4Constraints.begin());
		DKinFitParticle* locKinFitParticle = locP4Constraint->Get_DefinedParticle();
		if(locKinFitParticle != NULL)
			locKinFitParticle->Set_Momentum(locP4Constraint->Get_InitP3Guess());
	}

	//decaying: p3
	for(locParticleIterator = dKinFitParticles.begin(); locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		if(locKinFitParticle->Get_KinFitParticleType() == d_DecayingParticle)
			locKinFitParticle->Set_Momentum(dKinFitUtils->Calc_DecayingP4_ByPosition(locKinFitParticle, true).Vect());
	}

	//set vertex constraint flags: used if not accelerating
	for(locVertexIterator = locVertexConstraints.begin(); locVertexIterator != locVertexConstraints.end(); ++locVertexIterator)
	{
		set<DKinFitParticle*> locFullConstrainParticles = (*locVertexIterator)->Get_FullConstrainParticles();
		for(locParticleIterator = locFullConstrainParticles.begin(); locParticleIterator != locFullConstrainParticles.end(); ++locParticleIterator)
		{
			DKinFitParticle* locParticle = *locParticleIterator;
			TVector3 locMomentum = locParticle->Get_Momentum();

			unsigned short int locVertexConstraintFlag = 0;
			if(fabs(locMomentum.Pz()) > fabs(locMomentum.Px()))
				locVertexConstraintFlag = (fabs(locMomentum.Pz()) > fabs(locMomentum.Py())) ? 1 : 2;
			else
				locVertexConstraintFlag = (fabs(locMomentum.Px()) > fabs(locMomentum.Py())) ? 3 : 2;
			locParticle->Set_VertexConstraintFlag(locVertexConstraintFlag);
		}
	}
}

void DKinFitter::Set_MatrixSizes(void)
{
	//set matrix sizes
	dNumXi = 0; //num unknowns
	dNumEta = 0; //num measurables
	dNumF = 0; //num constraint eqs

	//Calculate dNumEta
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		bool locIsInP4Constraint = Get_IsInConstraint<DKinFitConstraint_P4>(locKinFitParticle);
		bool locIsInMassConstraint = Get_IsInConstraint<DKinFitConstraint_Mass>(locKinFitParticle);
		bool locIsInVertexConstraint = Get_IsInConstraint<DKinFitConstraint_Vertex>(locKinFitParticle);
		bool locIsIndirectlyInVertexConstraint = Get_IsIndirectlyInConstraint<DKinFitConstraint_Vertex>(locKinFitParticle);
		bool locChargedBFieldFlag = (locKinFitParticle->Get_Charge() != 0) && dKinFitUtils->Get_IsBFieldNearBeamline();

		if(dDebugLevel >= 20)
			cout << "PID, pointer, is in p4/mass/vert/indirectv, accel = " << locKinFitParticle->Get_PID() << ", " << locKinFitParticle << ", " << locIsInP4Constraint << ", " << locIsInMassConstraint << ", " << locIsInVertexConstraint << ", " << locIsIndirectlyInVertexConstraint << ", " << locChargedBFieldFlag << endl;
		if(!locKinFitParticle->Get_IsNeutralShowerFlag())
		{
			if(locIsInP4Constraint || locIsInMassConstraint || Get_IsVertexConstrained(locKinFitParticle) || locIsIndirectlyInVertexConstraint)
				dNumEta += 3; //p3
			if(Get_IsVertexConstrained(locKinFitParticle) || (locIsIndirectlyInVertexConstraint && locChargedBFieldFlag))
				dNumEta += 3; //v3 //directly (first condition) or indirectly AND accelerating (second condition)
		}
		else //neutral shower
		{
			if((locIsInP4Constraint || locIsInMassConstraint) && locIsInVertexConstraint)
				dNumEta += 4; //E + v3 (p4 fit needs p3, which is derived from v3 + kinfit vertex)
		}
		if(Get_IsTimeConstrained(locKinFitParticle))
			++dNumEta; //t //directly (first condition) or indirectly AND accelerating (second condition)
	}

	//Calculate dNumXi and dNumF
	set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
	for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(*locConstraintIterator);
		if(locKinFitConstraint_P4 != NULL)
		{
			dNumF += 4;
			if(locKinFitConstraint_P4->Get_DefinedParticle() != NULL)
				dNumXi += 3; //p3 //missing or open-ended decaying particles
			continue;
		}

		DKinFitConstraint_Mass* locKinFitConstraint_Mass = dynamic_cast<DKinFitConstraint_Mass*>(*locConstraintIterator);
		if(locKinFitConstraint_Mass != NULL)
		{
			dNumF += 1;
			continue;
		}

		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(*locConstraintIterator);
		if(locKinFitConstraint_Vertex != NULL)
		{
			dNumXi += 3; //v3
			dNumF += 2*locKinFitConstraint_Vertex->Get_FullConstrainParticles().size();
			continue;
		}

		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(*locConstraintIterator);
		if(locKinFitConstraint_Spacetime != NULL)
		{
			dNumXi += 4; //v3, t
			dNumF += 3*locKinFitConstraint_Spacetime->Get_FullConstrainParticles().size();
			dNumF += locKinFitConstraint_Spacetime->Get_OnlyConstrainTimeParticles().size(); //for each neutral shower
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

	//SETUP dY
	int locParamIndex = 0;
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		TVector3 locMomentum = locKinFitParticle->Get_Momentum();
		TVector3 locPosition = locKinFitParticle->Get_Position();

		bool locIsInP4Constraint = Get_IsInConstraint<DKinFitConstraint_P4>(locKinFitParticle);
		bool locIsInMassConstraint = Get_IsInConstraint<DKinFitConstraint_Mass>(locKinFitParticle);
		bool locIsInVertexConstraint = Get_IsInConstraint<DKinFitConstraint_Vertex>(locKinFitParticle);
		bool locIsIndirectlyInVertexConstraint = Get_IsIndirectlyInConstraint<DKinFitConstraint_Vertex>(locKinFitParticle);
		bool locChargedBFieldFlag = (locKinFitParticle->Get_Charge() != 0) && dKinFitUtils->Get_IsBFieldNearBeamline();

		if(!locKinFitParticle->Get_IsNeutralShowerFlag()) //non-neutral-shower
		{
			if(locIsInP4Constraint || locIsInMassConstraint || Get_IsVertexConstrained(locKinFitParticle) || locIsIndirectlyInVertexConstraint)
			{
				locKinFitParticle->Set_PxParamIndex(locParamIndex);
				dY(locParamIndex, 0) = locMomentum.Px();
				dY(locParamIndex + 1, 0) = locMomentum.Py();
				dY(locParamIndex + 2, 0) = locMomentum.Pz();
				locParamIndex += 3;
			}
			if(Get_IsVertexConstrained(locKinFitParticle) || (locIsIndirectlyInVertexConstraint && locChargedBFieldFlag))
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
			if((locIsInP4Constraint || locIsInMassConstraint) && locIsInVertexConstraint)
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

		if(Get_IsTimeConstrained(locKinFitParticle))
		{
			locKinFitParticle->Set_TParamIndex(locParamIndex);
			dY(locParamIndex, 0) = locKinFitParticle->Get_Time();
			++locParamIndex;
		}
	}

	//SETUP dEta
	dEta = dY; //use measurements as first guess

	//SETUP dXi (with initial guesses) and constraint equation indices
	locParamIndex = 0;
	int locConstraintIndex = 0;
	set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
	for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(*locConstraintIterator);
		if(locKinFitConstraint_P4 != NULL)
		{
			locKinFitConstraint_P4->Set_FIndex(locConstraintIndex);
			locConstraintIndex += 4; //p4

			DKinFitParticle* locDefinedParticle = locKinFitConstraint_P4->Get_DefinedParticle();
			if(locDefinedParticle == NULL)
				continue; //no missing particles

			//set initial p3 guess
			TVector3 locMomentum = locKinFitConstraint_P4->Get_InitP3Guess();
			dXi(locParamIndex, 0) = locMomentum.Px();
			dXi(locParamIndex + 1, 0) = locMomentum.Py();
			dXi(locParamIndex + 2, 0) = locMomentum.Pz();
			locDefinedParticle->Set_PxParamIndex(locParamIndex);
			locParamIndex += 3;
			continue;
		}

		DKinFitConstraint_Mass* locKinFitConstraint_Mass = dynamic_cast<DKinFitConstraint_Mass*>(*locConstraintIterator);
		if(locKinFitConstraint_Mass != NULL)
		{
			locKinFitConstraint_Mass->Set_FIndex(locConstraintIndex);
			locConstraintIndex += 1; //mass
			continue;
		}

		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(*locConstraintIterator);
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(*locConstraintIterator);
		if((locKinFitConstraint_Vertex != NULL) && (locKinFitConstraint_Spacetime == NULL))
		{
			TVector3 locPosition = locKinFitConstraint_Vertex->Get_InitVertexGuess();
			dXi(locParamIndex, 0) = locPosition.X();
			dXi(locParamIndex + 1, 0) = locPosition.Y();
			dXi(locParamIndex + 2, 0) = locPosition.Z();
			locKinFitConstraint_Vertex->Set_CommonVxParamIndex(locParamIndex);
			locParamIndex += 3;

			set<DKinFitParticle*> locFullConstrainParticles = locKinFitConstraint_Vertex->Get_FullConstrainParticles();
			locParticleIterator = locFullConstrainParticles.begin();
			for(; locParticleIterator != locFullConstrainParticles.end(); ++locParticleIterator)
			{
				locKinFitConstraint_Vertex->Set_FIndex(*locParticleIterator, locConstraintIndex);
				locConstraintIndex += 2;
			}
			continue;
		}

		if(locKinFitConstraint_Spacetime != NULL)
		{
			TVector3 locPosition = locKinFitConstraint_Spacetime->Get_InitVertexGuess();
			dXi(locParamIndex, 0) = locPosition.X();
			dXi(locParamIndex + 1, 0) = locPosition.Y();
			dXi(locParamIndex + 2, 0) = locPosition.Z();
			dXi(locParamIndex + 3, 0) = locKinFitConstraint_Spacetime->Get_InitTimeGuess();
			locKinFitConstraint_Spacetime->Set_CommonVxParamIndex(locParamIndex);
			locKinFitConstraint_Spacetime->Set_CommonTParamIndex(locParamIndex + 3);
			locParamIndex += 4;

			set<DKinFitParticle*> locFullConstrainParticles = locKinFitConstraint_Spacetime->Get_FullConstrainParticles();
			locParticleIterator = locFullConstrainParticles.begin();
			for(; locParticleIterator != locFullConstrainParticles.end(); ++locParticleIterator)
			{
				locKinFitConstraint_Spacetime->Set_FIndex(*locParticleIterator, locConstraintIndex);
				locConstraintIndex += 3;
			}

			set<DKinFitParticle*> locOnlyConstrainTimeParticles = locKinFitConstraint_Spacetime->dOnlyConstrainTimeParticles;
			locParticleIterator = locOnlyConstrainTimeParticles.begin();
			for(; locParticleIterator != locOnlyConstrainTimeParticles.end(); ++locParticleIterator)
			{
				locKinFitConstraint_Spacetime->Set_FIndex(*locParticleIterator, locConstraintIndex);
				++locConstraintIndex; 
			}

			continue;
		}
	}

	//SETUP dVY
	locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue; //uncertainties of target particle momentum are zero

		TVector3 locMomentum = locKinFitParticle->Get_Momentum();
		TVector3 locPosition = locKinFitParticle->Get_Position();
		const TMatrixDSym& locCovarianceMatrix = *(locKinFitParticle->Get_CovarianceMatrix());

		int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
		int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
		int locTParamIndex = locKinFitParticle->Get_TParamIndex();
		int locEParamIndex = locKinFitParticle->Get_EParamIndex();

		int locCovMatrixEParamIndex = locKinFitParticle->Get_CovMatrixEParamIndex();
		int locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
		int locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
		int locCovMatrixTParamIndex = locKinFitParticle->Get_CovMatrixTParamIndex();

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

	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dEta: " << endl;
		Print_Matrix(dEta);
		cout << "DKinFitter: dVY: " << endl;
		Print_Matrix(dVY);
		cout << "DKinFitter: dXi: " << endl;
		Print_Matrix(dXi);
		for(locParticleIterator = dKinFitParticles.begin(); locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
			(*locParticleIterator)->Print_ParticleParams();
	}
}

/********************************************************************* PERFORM FIT *********************************************************************/

bool DKinFitter::Fit_Reaction(void)
{
	//Initialize
	Prepare_ConstraintsAndParticles();
	if(!dKinFitUtils->Validate_Constraints(dKinFitConstraints))
	{
		dKinFitStatus = d_KinFitFailedSetup;
		return false;
	}
	Set_MatrixSizes();
	Resize_Matrices();
	Fill_InputMatrices();

	//Do the fit
	if(!Iterate())
		return false;

	// Calculate final covariance matrices
	if(dNumXi > 0)
		*dVXi = dU;
	Calc_dVdEta();

	// Calculate fit NDF, Confidence level
	dNDF = dNumF - dNumXi;
	dConfidenceLevel = TMath::Prob(dChiSq, dNDF);

	// Calculate pulls
	dEpsilon = dY - dEta;
	Calc_Pulls();

	// Set final track info
	Set_FinalTrackInfo(); //Update_ParticleParams() already called: at the end of the last iteration 

	if(dDebugLevel > 5)
		cout << "DKinFitter: Final dChiSq, dNDF, dConfidenceLevel = " << dChiSq << ", " << dNDF << ", " << dConfidenceLevel << endl;

	dKinFitStatus = d_KinFitSuccessful;
	return true;
}

bool DKinFitter::Iterate(void)
{
	dChiSq = 9.9E99;
	double locPreviousChiSq = 0.0;
	unsigned int locNumIterations = 0;
	TMatrixD locR(dNumF, 1);
	while((fabs(dChiSq - locPreviousChiSq) > dConvergenceChiSqDiff) || (dChiSq < 0.0))
	{
		if(locNumIterations >= dMaxNumIterations)
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
			dKinFitStatus = d_KinFitFailedInversion;
			return false; // matrix is not invertible
		}

		if(dNumXi > 0)
		{
			if(!Calc_dU())
			{
				dKinFitStatus = d_KinFitFailedInversion;
				return false; // matrix is not invertible
			}

			TMatrixD locDeltaXi(dNumXi, 1);
			locDeltaXi = -1.0*dU*dF_dXi_T*dS_Inverse*locR;

			dXi += locDeltaXi;
			if(dDebugLevel > 20)
			{
				cout << "DKinFitter: locDeltaXi: " << endl;
				Print_Matrix(locDeltaXi);
				cout << "DKinFitter: dXi: " << endl;
				Print_Matrix(dXi);
			}

			dLambda = dS_Inverse*(locR + dF_dXi*locDeltaXi);
		}
		else
			dLambda = dS_Inverse*locR;

		dLambda_T.Transpose(dLambda);
		dEta = dY - dVY*dF_dEta_T*dLambda;

		TMatrixDSym locTempMatrix = dS; //similarity (below) destroys the matrix: use a temp to preserve dS
		dChiSq = (locTempMatrix.SimilarityT(dLambda) + 2.0*dLambda_T*dF)(0, 0);

		Update_ParticleParams(); //input eta & xi info into particle objects

		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dLambda: " << endl;
			Print_Matrix(dLambda);
			cout << "DKinFitter: dEta: " << endl;
			Print_Matrix(dEta);
			cout << "DKinFitter: dChiSq = " << dChiSq << endl;
			set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
			for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
				(*locParticleIterator)->Print_ParticleParams();
		}

		++locNumIterations;
	}

	return true;
}

/****************************************************************** CALCULATE MATRICES *****************************************************************/

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

void DKinFitter::Calc_dVdEta(void)
{
	TMatrixDSym locG = dS_Inverse;
	locG.SimilarityT(dF_dEta);
	if(dNumXi == 0)
	{
		*dVEta = dVY - locG.Similarity(dVY); //destroys locG, but it's not needed anymore
		*dV = *dVEta;
		if(dDebugLevel > 20)
		{
			cout << "DKinFitter: dVEta: " << endl;
			Print_Matrix(*dVEta);
			cout << "DKinFitter: dV: " << endl;
			Print_Matrix(*dV);
		}
		return;
	}

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

	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dVEta: " << endl;
		Print_Matrix(*dVEta);
		cout << "DKinFitter: dV: " << endl;
		Print_Matrix(*dV);
	}
}

void DKinFitter::Calc_dF(void)
{
	dF.Zero();
	dF_dXi.Zero();
	dF_dEta.Zero();

	size_t locFIndex = 0;
	set<DKinFitConstraint*>::iterator locConstraintIterator = dKinFitConstraints.begin();
	for(; locConstraintIterator != dKinFitConstraints.end(); ++locConstraintIterator)
	{
		DKinFitConstraint_P4* locKinFitConstraint_P4 = dynamic_cast<DKinFitConstraint_P4*>(*locConstraintIterator);
		if(locKinFitConstraint_P4 != NULL)
		{
			int locFIndex = locKinFitConstraint_P4->Get_FIndex();
			if(dDebugLevel > 10)
				cout << "DKinFitter: P4 Constraint, F index = " << locFIndex << endl;

			//initial state
			set<DKinFitParticle*> locInitialParticles = locKinFitConstraint_P4->Get_InitialParticles();
			set<DKinFitParticle*>::iterator locParticleIterator = locInitialParticles.begin();
			for(; locParticleIterator != locInitialParticles.end(); ++locParticleIterator)
				Calc_dF_P4(locFIndex, *locParticleIterator, 1.0);

			//final state
			set<DKinFitParticle*> locFinalParticles = locKinFitConstraint_P4->Get_FinalParticles();
			for(locParticleIterator = locFinalParticles.begin(); locParticleIterator != locFinalParticles.end(); ++locParticleIterator)
				Calc_dF_P4(locFIndex, *locParticleIterator, -1.0);

			continue;
		}

		DKinFitConstraint_Mass* locKinFitConstraint_Mass = dynamic_cast<DKinFitConstraint_Mass*>(*locConstraintIterator);
		if(locKinFitConstraint_Mass != NULL)
		{
			int locFIndex = locKinFitConstraint_Mass->Get_FIndex();
			if(dDebugLevel > 10)
				cout << "DKinFitter: Mass Constraint, F index = " << locFIndex << endl;

			DKinFitParticle* locDecayingParticle = locKinFitConstraint_Mass->Get_DecayingParticle();
			double locTargetedMass = locDecayingParticle->Get_Mass();

			TLorentzVector locXP4 = dKinFitUtils->Calc_DecayingP4_ByP3Derived(locDecayingParticle, true);
			if(dDebugLevel > 30)
				cout << "Final decaying pxyzE is: " << locXP4.Px() << ", " << locXP4.Py() << ", " << locXP4.Pz() << ", " << locXP4.E() << endl;
			dF(locFIndex, 0) += locXP4.M2() - locTargetedMass*locTargetedMass;

			Calc_dF_MassDerivs(locFIndex, locDecayingParticle, locXP4, 1.0, true);
			continue;
		}

		DKinFitConstraint_Vertex* locKinFitConstraint_Vertex = dynamic_cast<DKinFitConstraint_Vertex*>(*locConstraintIterator);
		DKinFitConstraint_Spacetime* locKinFitConstraint_Spacetime = dynamic_cast<DKinFitConstraint_Spacetime*>(*locConstraintIterator);
		if((locKinFitConstraint_Vertex != NULL) && (locKinFitConstraint_Spacetime == NULL))
		{
			set<DKinFitParticle*> locFullConstrainParticles = locKinFitConstraint_Vertex->Get_FullConstrainParticles();
			set<DKinFitParticle*>::iterator locParticleIterator = locFullConstrainParticles.begin();
			for(; locParticleIterator != locFullConstrainParticles.end(); ++locParticleIterator)
			{
				DKinFitParticle* locKinFitParticle = *locParticleIterator;
				DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

				bool locIsDecayingFlag = (locKinFitParticleType == d_DecayingParticle);
				locFIndex = locKinFitConstraint_Vertex->Get_FIndex(locKinFitParticle);
				if(dDebugLevel > 10)
					cout << "DKinFitter: F index, locIsDecayingFlag = " << locFIndex << ", " << locIsDecayingFlag << endl;

				Calc_dF_Vertex(locFIndex, locKinFitParticle, NULL, 1.0);
			}

			continue;
		}

		if(locKinFitConstraint_Spacetime != NULL)
		{
			//HANDLE HERE WHEN READY!
		}
	}

	dF_dEta_T.Transpose(dF_dEta);
	dF_dXi_T.Transpose(dF_dXi);
}

/**************************************************************** CALCULATE P4 MATRICES ****************************************************************/

void DKinFitter::Calc_dF_P4(int locFIndex, const DKinFitParticle* locKinFitParticle, double locStateSignMultiplier)
{
	//E, px, py, pz
	int locCharge = locKinFitParticle->Get_Charge();
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locBField = dKinFitUtils->Get_IsBFieldNearBeamline() ? dKinFitUtils->Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
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
	bool locChargedBFieldFlag = (locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline();
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locKinFitParticleType != d_DecayingParticle) || (locPxParamIndex >= 0))
	{
		//not an enclosed decaying particle. for decaying particles, will instead get p4 from the deriving particles
		dF(locFIndex, 0) += locStateSignMultiplier*locP4.E();
		dF(locFIndex + 1, 0) += locStateSignMultiplier*locP4.Px();
		dF(locFIndex + 2, 0) += locStateSignMultiplier*locP4.Py();
		dF(locFIndex + 3, 0) += locStateSignMultiplier*locP4.Pz();
	}

	if(dDebugLevel > 30)
		cout << "PID, sign, pxyzE = " << locKinFitParticle->Get_PID() << ", " << locStateSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(locCommonVertexFitFlag && locChargedBFieldFlag && (locKinFitParticleType != d_TargetParticle) && (locKinFitParticleType != d_MissingParticle))
	{
		//fitting vertex of charged track in magnetic field: momentum changes as function of vertex
		//decaying particles: p4 not directly used, deriving particles are: so must propagate if charged
		//if initial particle is a "detected" particle (actually a decaying particle treated as detected): still propagate vertex (assume p3/v3 defined at production vertex)

		TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
		if(dDebugLevel > 30)
			cout << "propagate pxyz by: " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.X() << ", " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.Y() << ", " << -1.0*locStateSignMultiplier*locA*locDeltaXCrossH.Z() << endl;
		dF(locFIndex + 1, 0) -= locStateSignMultiplier*locA*locDeltaXCrossH.X();
		dF(locFIndex + 2, 0) -= locStateSignMultiplier*locA*locDeltaXCrossH.Y();
		dF(locFIndex + 3, 0) -= locStateSignMultiplier*locA*locDeltaXCrossH.Z();
	}

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 1; PID = " << locKinFitParticle->Get_PID() << endl;
		return; //target params are fixed: no partial derivatives
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 3; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dEta(locFIndex, locEParamIndex) = locStateSignMultiplier;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dF_dEta(locFIndex + 1, locEParamIndex) = locStateSignMultiplier*locEOverPSq*locP4.Px();
		dF_dEta(locFIndex + 2, locEParamIndex) = locStateSignMultiplier*locEOverPSq*locP4.Py();
		dF_dEta(locFIndex + 3, locEParamIndex) = locStateSignMultiplier*locEOverPSq*locP4.Pz();

		TVector3 locDeltaXOverMagDeltaXSq = locDeltaX*(1.0/locDeltaX.Mag2());

		dF_dEta(locFIndex + 1, locVxParamIndex) = locStateSignMultiplier*locP4.Px()*(locDeltaXOverMagDeltaXSq.X() - 1.0/locDeltaX.X());
		dF_dEta(locFIndex + 2, locVxParamIndex + 1) = locStateSignMultiplier*locP4.Py()*(locDeltaXOverMagDeltaXSq.Y() - 1.0/locDeltaX.Y());
		dF_dEta(locFIndex + 3, locVxParamIndex + 2) = locStateSignMultiplier*locP4.Pz()*(locDeltaXOverMagDeltaXSq.Z() - 1.0/locDeltaX.Z());

		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locStateSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Y();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locStateSignMultiplier*locP4.Px()*locDeltaXOverMagDeltaXSq.Z();

		dF_dEta(locFIndex + 2, locVxParamIndex) = locStateSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.X();
		dF_dEta(locFIndex + 2, locVxParamIndex + 2) = locStateSignMultiplier*locP4.Py()*locDeltaXOverMagDeltaXSq.Z();

		dF_dEta(locFIndex + 3, locVxParamIndex) = locStateSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.X();
		dF_dEta(locFIndex + 3, locVxParamIndex + 1) = locStateSignMultiplier*locP4.Pz()*locDeltaXOverMagDeltaXSq.Y();

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
		//missing particle or open-ended decaying particle: p3 is unknown
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 4; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dXi(locFIndex, locPxParamIndex) = locStateSignMultiplier*locP4.Px()/locP4.E();
		dF_dXi(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locP4.Py()/locP4.E();
		dF_dXi(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locP4.Pz()/locP4.E();

		dF_dXi(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier;
		dF_dXi(locFIndex + 2, locPxParamIndex + 1) = locStateSignMultiplier;
		dF_dXi(locFIndex + 3, locPxParamIndex + 2) = locStateSignMultiplier;
	}
	else if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 5; PID = " << locKinFitParticle->Get_PID() << endl;
		if(locChargedBFieldFlag && locCommonVertexFitFlag)
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_P4() Section 5a" << endl;
			//vertex factors
				//locVertexSignMultiplier is needed because the signs for Vx & CommonVx might be switched
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locStateSignMultiplier*locVertexSignMultiplier*locA*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += -1.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locH.Y();

			dF_dXi(locFIndex + 2, locVxParamIndex) += -1.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locH.Z();
			dF_dXi(locFIndex + 2, locVxParamIndex + 2) += locStateSignMultiplier*locVertexSignMultiplier*locA*locH.X();

			dF_dXi(locFIndex + 3, locVxParamIndex) += locStateSignMultiplier*locVertexSignMultiplier*locA*locH.Y();
			dF_dXi(locFIndex + 3, locVxParamIndex + 1) += -1.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locH.X();

			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);

			dF_dXi(locFIndex + 2, locCommonVxParamIndex) -= dF_dXi(locFIndex + 2, locVxParamIndex);
			dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 2, locVxParamIndex + 2);

			dF_dXi(locFIndex + 3, locCommonVxParamIndex) -= dF_dXi(locFIndex + 3, locVxParamIndex);
			dF_dXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 3, locVxParamIndex + 1);
		}

		//replace the decaying particle with the particles it's momentum is derived from
		//initial state
		set<DKinFitParticle*> locFromInitialState = locKinFitParticle->Get_FromInitialState();
		set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
		for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			Calc_dF_P4(locFIndex, *locParticleIterator, locStateSignMultiplier); //decaying particle multiplier * 1.0
		}

		//final state
		set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		//If defined by invariant mass: add p4s of final state particles
		//If defined by missing mass: add p4s of init state, subtract final state
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			Calc_dF_P4(locFIndex, *locParticleIterator, locNextStateSignMultiplier);
		}
	}
	else if(locChargedBFieldFlag && locCommonVertexFitFlag && (locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
	{
		//detected charged particle in b-field (can be beam particle) & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 2; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = locStateSignMultiplier*locP4.Px()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locP4.Py()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locP4.Pz()/locP4.E();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier;
		dF_dEta(locFIndex + 2, locPxParamIndex + 1) = locStateSignMultiplier;
		dF_dEta(locFIndex + 3, locPxParamIndex + 2) = locStateSignMultiplier;

		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locStateSignMultiplier*locA*locH.Z();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locStateSignMultiplier*locA*locH.Y();

		dF_dEta(locFIndex + 2, locVxParamIndex) = -1.0*locStateSignMultiplier*locA*locH.Z();
		dF_dEta(locFIndex + 2, locVxParamIndex + 2) = locStateSignMultiplier*locA*locH.X();

		dF_dEta(locFIndex + 3, locVxParamIndex) = locStateSignMultiplier*locA*locH.Y();
		dF_dEta(locFIndex + 3, locVxParamIndex + 1) = -1.0*locStateSignMultiplier*locA*locH.X();

		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);

		dF_dXi(locFIndex + 2, locCommonVxParamIndex) -= dF_dEta(locFIndex + 2, locVxParamIndex);
		dF_dXi(locFIndex + 2, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 2, locVxParamIndex + 2);

		dF_dXi(locFIndex + 3, locCommonVxParamIndex) -= dF_dEta(locFIndex + 3, locVxParamIndex);
		dF_dXi(locFIndex + 3, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 3, locVxParamIndex + 1);
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_P4() Section 6; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = locStateSignMultiplier*locP4.Px()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locP4.Py()/locP4.E();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locP4.Pz()/locP4.E();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier;
		dF_dEta(locFIndex + 2, locPxParamIndex + 1) = locStateSignMultiplier;
		dF_dEta(locFIndex + 3, locPxParamIndex + 2) = locStateSignMultiplier;
	}
}

/*************************************************************** CALCULATE MASS MATRICES ***************************************************************/

void DKinFitter::Calc_dF_MassDerivs(size_t locFIndex, const DKinFitParticle* locKinFitParticle, TLorentzVector locXP4, double locStateSignMultiplier, bool locIsConstrainedParticle)
{
	//E, px, py, pz
	int locCharge = locKinFitParticle->Get_Charge();
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	TLorentzVector locP4 = locKinFitParticle->Get_P4();
	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locBField = dKinFitUtils->Get_IsBFieldNearBeamline() ? dKinFitUtils->Get_BField(locPosition) : TVector3(0.0, 0.0, 0.0);
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
		//Thus we need a factor of -1 for the Xi-
	bool locNeedP4AtProductionVertex = locKinFitParticle->Get_FromInitialState().empty(); //true if defined by decay products; else by missing mass
	double locVertexSignMultiplier = (locNeedP4AtProductionVertex == locKinFitParticle->Get_VertexP4AtProductionVertex()) ? -1.0 : 1.0;
	TVector3 locDeltaX = locVertexSignMultiplier*(locCommonVertex - locPosition); //vector points in the OPPOSITE direction of the momentum

	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	bool locCommonVertexFitFlag = locKinFitParticle->Get_FitCommonVertexFlag();
	bool locChargedBFieldFlag = (locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline();
	bool locNeutralShowerFlag = locKinFitParticle->Get_IsNeutralShowerFlag();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locEParamIndex = locKinFitParticle->Get_EParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	TVector3 locPXCrossH = locXP4.Vect().Cross(locH);

	if(dDebugLevel > 30)
		cout << "PID, sign, pxyzE = " << locKinFitParticle->Get_PID() << ", " << locStateSignMultiplier << ", " << locP4.Px() << ", " << locP4.Py() << ", " << locP4.Pz() << ", " << locP4.E() << endl;

	if(locKinFitParticleType == d_TargetParticle)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 1; PID = " << locKinFitParticle->Get_PID() << endl;
		return; //target params are fixed: no partial derivatives
	}
	else if(locNeutralShowerFlag)
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 3; PID = " << locKinFitParticle->Get_PID() << endl;

		double locEOverPSq = locP4.E()/locP4.Vect().Mag2();
		dF_dEta(locFIndex, locEParamIndex) = 2.0*locStateSignMultiplier*(locXP4.E() - locEOverPSq*locXP4.Vect().Dot(locP4.Vect()));

		double locDeltaXDotPXOverMagDeltaXSq = locDeltaX.Dot(locXP4.Vect())/(locDeltaX.Mag2());
		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locStateSignMultiplier*locP4.Px()*(locXP4.Px()/locDeltaX.X() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locStateSignMultiplier*locP4.Py()*(locXP4.Py()/locDeltaX.Y() - locDeltaXDotPXOverMagDeltaXSq);
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locStateSignMultiplier*locP4.Pz()*(locXP4.Pz()/locDeltaX.Z() - locDeltaXDotPXOverMagDeltaXSq);

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else if(locKinFitParticleType == d_MissingParticle)
	{
		//missing or open-ended-decaying particle: p3 is unknown (not derivable)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 4; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dXi(locFIndex, locPxParamIndex) = 2.0*locStateSignMultiplier*(locP4.Px()*locXP4.E()/locP4.E() - locXP4.Px());
		dF_dXi(locFIndex, locPxParamIndex + 1) = 2.0*locStateSignMultiplier*(locP4.Py()*locXP4.E()/locP4.E() - locXP4.Py());
		dF_dXi(locFIndex, locPxParamIndex + 2) = 2.0*locStateSignMultiplier*(locP4.Pz()*locXP4.E()/locP4.E() - locXP4.Pz());
	}
	else if(locKinFitParticleType == d_DecayingParticle)
	{
		//enclosed decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 5; PID = " << locKinFitParticle->Get_PID() << endl;
		if(locChargedBFieldFlag && locCommonVertexFitFlag && !locIsConstrainedParticle)
		{
			//if locKinFitSubConstraint_P4 is NULL is first particle; don't include: replace it
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_MassDerivs() Section 5a" << endl;

			dF_dXi(locFIndex, locVxParamIndex) += 2.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locPXCrossH.X();
			dF_dXi(locFIndex, locVxParamIndex + 1) += 2.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locPXCrossH.Y();
			dF_dXi(locFIndex, locVxParamIndex + 2) += 2.0*locStateSignMultiplier*locVertexSignMultiplier*locA*locPXCrossH.Z();

			dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);
		}

		//replace the decaying particle with the particles it's momentum is derived from
		//initial state
		set<DKinFitParticle*> locFromInitialState = locKinFitParticle->Get_FromInitialState();
		set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
		for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with init-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			Calc_dF_MassDerivs(locFIndex, *locParticleIterator, locXP4, locStateSignMultiplier, false); //decaying particle multiplier * 1.0
		}

		//final state
		set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
		bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
		double locNextStateSignMultiplier = locStateSignMultiplier;
		//If defined by invariant mass: add p4s of final state particles
		//If defined by missing mass: add p4s of init state, subtract final state
		if(!locDefinedByInvariantMassFlag)
			locNextStateSignMultiplier *= -1.0;
		for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
		{
			if(dDebugLevel > 30)
				cout << "decaying, partially replace with final-state PID = " << (*locParticleIterator)->Get_PID() << endl;
			Calc_dF_MassDerivs(locFIndex, *locParticleIterator, locXP4, locNextStateSignMultiplier, false);
		}
	}
	else if(locChargedBFieldFlag && locCommonVertexFitFlag && (locKinFitParticleType != d_DecayingParticle) && (locKinFitParticleType != d_MissingParticle))
	{
		//detected charged particle in b-field (can be beam particle) & in vertex fit
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 2; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*locStateSignMultiplier*(locP4.Px()*locXP4.E()/locP4.E() - locXP4.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*locStateSignMultiplier*(locP4.Py()*locXP4.E()/locP4.E() - locXP4.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*locStateSignMultiplier*(locP4.Pz()*locXP4.E()/locP4.E() - locXP4.Pz());

		dF_dEta(locFIndex, locVxParamIndex) = 2.0*locStateSignMultiplier*locA*locPXCrossH.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = 2.0*locStateSignMultiplier*locA*locPXCrossH.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = 2.0*locStateSignMultiplier*locA*locPXCrossH.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);
	}
	else
	{
		// either no common vertex constraint, charged and detected but b-field = 0, or neutral particle with pre-ordained vertex (e.g. beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_MassDerivs() Section 6; PID = " << locKinFitParticle->Get_PID() << endl;

		dF_dEta(locFIndex, locPxParamIndex) = 2.0*locStateSignMultiplier*(locP4.Px()*locXP4.E()/locP4.E() - locXP4.Px());
		dF_dEta(locFIndex, locPxParamIndex + 1) = 2.0*locStateSignMultiplier*(locP4.Py()*locXP4.E()/locP4.E() - locXP4.Py());
		dF_dEta(locFIndex, locPxParamIndex + 2) = 2.0*locStateSignMultiplier*(locP4.Pz()*locXP4.E()/locP4.E() - locXP4.Pz());
	}
}

/************************************************************** CALCULATE VERTEX MATRICES **************************************************************/

void DKinFitter::Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier)
{
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

	if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 1; PID = " << locKinFitParticle->Get_PID() << endl;
		return; //no partial derivatives
	}

	if(locKinFitParticle_DecayingSource == NULL) //this particle is directly at the vertex
		Calc_dF_Vertex_NotDecaying(locFIndex, locKinFitParticle);
	else //this particle is used to define the momentum of a decaying particle
	{
		int locCharge = locKinFitParticle_DecayingSource->Get_Charge();
		if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline())
			Calc_dF_Vertex_Decaying_Accel(locFIndex, locKinFitParticle, locKinFitParticle_DecayingSource, locStateSignMultiplier);
		else
			Calc_dF_Vertex_Decaying_NonAccel(locFIndex, locKinFitParticle, locKinFitParticle_DecayingSource, locStateSignMultiplier);
	}

	if(locKinFitParticleType != d_DecayingParticle)
		return;

	// find the p4 constraint where the decaying particle momentum is defined //this is the point where the position is defined as well
		//why this constraint? because if you choose the other one, the momentum has to be propagated across a distance, requiring a different (& MUCH more complicated) partial derivative

	//replace the decaying particle with the particles it's momentum is derived from
	//initial state
	set<DKinFitParticle*> locFromInitialState = locKinFitParticle->Get_FromInitialState();
	set<DKinFitParticle*>::iterator locParticleIterator = locFromInitialState.begin();
	for(; locParticleIterator != locFromInitialState.end(); ++locParticleIterator)
	{
		if(dDebugLevel > 30)
			cout << "decaying, partially replace with init-state PID = " << (*locParticleIterator)->Get_PID() << endl;
		Calc_dF_Vertex(locFIndex, *locParticleIterator, locKinFitParticle, locStateSignMultiplier);
	}

	//final state
	set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
	bool locDefinedByInvariantMassFlag = locFromInitialState.empty();
	double locNextStateSignMultiplier = locStateSignMultiplier;
	//If defined by invariant mass: add p4s of final state particles
	//If defined by missing mass: add p4s of init state, subtract final state
	if(!locDefinedByInvariantMassFlag)
		locNextStateSignMultiplier *= -1.0;
	for(locParticleIterator = locFromFinalState.begin(); locParticleIterator != locFromFinalState.end(); ++locParticleIterator)
	{
		if(dDebugLevel > 30)
			cout << "decaying, partially replace with final-state PID = " << (*locParticleIterator)->Get_PID() << endl;
		Calc_dF_Vertex(locFIndex, *locParticleIterator, locKinFitParticle, locNextStateSignMultiplier);
	}
}

void DKinFitter::Calc_dF_Vertex_NotDecaying(size_t locFIndex, const DKinFitParticle* locKinFitParticle)
{
	//The input particle is directly at the vertex being constrained
		//The input particle ITSELF may be decaying
	DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
	int locCharge = locKinFitParticle->Get_Charge();

	TVector3 locPosition = locKinFitParticle->Get_Position();
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex();
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locMomentum = locKinFitParticle->Get_Momentum();

	int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
	int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
	int locCommonVxParamIndex = locKinFitParticle->Get_CommonVxParamIndex();

	if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 2; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
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
	else if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle))
	{
		//constraining this decaying charged particle in b-field //one-time contributions from decaying particle, does not include the particles it is replaced by (elsewhere)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 3; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle, locJ, locQ, locM, locD);

		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locPCrossH = locMomentum.Cross(locH);
		TVector3 locPCrossDeltaX = locMomentum.Cross(locDeltaX);
		double locDeltaXDotH = locDeltaX.Dot(locH);
		double locPDotH = locMomentum.Dot(locH);

		dF(locFIndex, 0) = locPCrossDeltaX.Dot(locH) - 0.5*locA*(locDeltaX.Mag2() - locDeltaXDotH*locDeltaXDotH);
		dF(locFIndex + 1, 0) = locDeltaXDotH - (locPDotH/locA)*asin(locJ);

		//true if the object's p3, & x4 are defined at its production vertex (& common x4 is at decay vertex). 
		bool locVertexP4AtProductionVertexFlag = locKinFitParticle->Get_VertexP4AtProductionVertex();

		//true if defined by decay products
		bool locP4DefinedByInvariantMassFlag = locKinFitParticle->Get_FromInitialState().empty();

		//Tricky when: P4-define step = common-vertex step
			//In other words: P4 calculated by particles at a different vertex where the position is defined
			//Main case: The Sigma+ in:   g, p -> K0, Sigma+     K0 -> pi+, pi-    Sigma+ -> p, pi0
		bool locP4DerivedAtCommonVertexFlag = (locP4DefinedByInvariantMassFlag == locVertexP4AtProductionVertexFlag);

		if(locP4DerivedAtCommonVertexFlag) //Tricky case
		{
			TVector3 locR = Calc_VertexParams_P4DerivedAtCommonVertex(locKinFitParticle);

			dF_dXi(locFIndex, locVxParamIndex) += locPCrossH.X();
			dF_dXi(locFIndex, locVxParamIndex + 1) += locPCrossH.Y();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locPCrossH.Z();

			dF_dXi(locFIndex + 1, locVxParamIndex) += locR.X();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locR.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += locR.Z();
		}
		else
		{
			dF_dXi(locFIndex, locVxParamIndex) += locPCrossH.X() + locA*locM.X();
			dF_dXi(locFIndex, locVxParamIndex + 1) += locPCrossH.Y() + locA*locM.Y();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locPCrossH.Z() + locA*locM.Z();

			dF_dXi(locFIndex + 1, locVxParamIndex) += locD.X();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locD.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += locD.Z();
		}

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dXi(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticleType == d_DecayingParticle) //non-accel decaying particle
	{
		//constraining this decaying non-accel particle //one-time contributions from decaying particle, does not include the particles it is replaced by (elsewhere)

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex + 1) += locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) += -1.0*locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) += -1.0*locMomentum.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) += -1.0*dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 2);
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex + 1) += locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) += -1.0*locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) += -1.0*dF_dXi(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		}
		else //2 & 3 //px is largest
		{
			dF(locFIndex, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 4c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex) += -1.0*locMomentum.Z();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locMomentum.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locMomentum.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex) += -1.0*dF_dXi(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dXi(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) += -1.0*dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		}
	}
	else //neutral detected or beam particle, or no magnetic field
	{
		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) += -1.0*dF_dEta(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 2);
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			dF(locFIndex, 0) = locMomentum.Y()*locDeltaX.Z() - locMomentum.Z()*locDeltaX.Y();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex + 1) = locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex + 1) += -1.0*dF_dEta(locFIndex, locVxParamIndex + 1);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		}
		else //2 & 3 //px is largest
		{
			dF(locFIndex, 0) = locMomentum.Z()*locDeltaX.X() - locMomentum.X()*locDeltaX.Z();
			dF(locFIndex + 1, 0) = locMomentum.X()*locDeltaX.Y() - locMomentum.Y()*locDeltaX.X();

			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 5c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locDeltaX.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locDeltaX.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locMomentum.Z();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locMomentum.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locMomentum.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locMomentum.X();

			dF_dXi(locFIndex, locCommonVxParamIndex) += -1.0*dF_dEta(locFIndex, locVxParamIndex);
			dF_dXi(locFIndex, locCommonVxParamIndex + 2) += -1.0*dF_dEta(locFIndex, locVxParamIndex + 2);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex);
			dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) += -1.0*dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		}
	}
}

void DKinFitter::Calc_dF_Vertex_Decaying_Accel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier)
{
	//The input particle is used to (perhaps indirectly) define a decaying particle at the vertex being constrained
		//The decaying particle that this particle directly defines is accelerating: charged and in a magnetic field
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

	if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 6; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locBField_DecayingSource = dKinFitUtils->Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();
		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		TVector3 locDeltaXCrossH_DecayingSource_CrossH = locDeltaXCrossH_DecayingSource.Cross(locH);
		TVector3 locQCrossH = locQ.Cross(locH);

		dF_dEta(locFIndex, locPxParamIndex) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dEta(locFIndex, locVxParamIndex) = -1.0*locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.X();
		dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Y();
		dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Z();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier*locQ.X();
		dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locStateSignMultiplier*locQ.Y();
		dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locStateSignMultiplier*locQ.Z();

		dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locStateSignMultiplier*locA*locQCrossH.X();
		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locStateSignMultiplier*locA*locQCrossH.Y();
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locStateSignMultiplier*locA*locQCrossH.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag() && (locKinFitParticle->Get_Mass() < 0.00001))
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 7; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locBField_DecayingSource = dKinFitUtils->Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		double locXTerm = locDeltaXCrossH_DecayingSource.Dot(locDeltaX)/(locDeltaX.Mag2());
		double locQTerm = locQ.Dot(locDeltaX)/(locDeltaX.Mag2());

		dF_dEta(locFIndex, locEParamIndex) = locStateSignMultiplier*locEnergy*(locDeltaXCrossH_DecayingSource.Dot(locMomentum))/(locMomentum.Mag2());
		dF_dEta(locFIndex, locVxParamIndex) = -1.0*locStateSignMultiplier*locMomentum.Px()*(locDeltaXCrossH_DecayingSource.X()/(locDeltaX.X()) - locXTerm);
		dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locStateSignMultiplier*locMomentum.Py()*(locDeltaXCrossH_DecayingSource.Y()/(locDeltaX.Y()) - locXTerm);
		dF_dEta(locFIndex, locVxParamIndex + 2) = -1.0*locStateSignMultiplier*locMomentum.Pz()*(locDeltaXCrossH_DecayingSource.Z()/(locDeltaX.Z()) - locXTerm);

		dF_dEta(locFIndex + 1, locEParamIndex) = locStateSignMultiplier*locEnergy*(locMomentum.Dot(locQ))/(locMomentum.Mag2());
		dF_dEta(locFIndex + 1, locVxParamIndex) = -1.0*locStateSignMultiplier*locMomentum.Px()*(locQ.X()/(locDeltaX.X()) - locQTerm);
		dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locStateSignMultiplier*locMomentum.Py()*(locQ.Y()/(locDeltaX.Y()) - locQTerm);
		dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locStateSignMultiplier*locMomentum.Pz()*(locQ.Z()/(locDeltaX.Z()) - locQTerm);

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag() && (locKinFitParticle->Get_Mass() >= 0.00001))
	{
		return;
	}
	else if((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle))
	{
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 8; PID = " << locKinFitParticle->Get_PID() << endl;

		//detected/beam, non-accelerating particle
		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locBField_DecayingSource = dKinFitUtils->Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);

		dF_dEta(locFIndex, locPxParamIndex) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dEta(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier*locQ.X();
		dF_dEta(locFIndex + 1, locPxParamIndex + 1) = locStateSignMultiplier*locQ.Y();
		dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locStateSignMultiplier*locQ.Z();
	}
	else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locKinFitParticle->Get_PxParamIndex() >= 0)))
	{
		//missing, or open-ended decaying particle
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 9; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locBField_DecayingSource = dKinFitUtils->Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);

		dF_dXi(locFIndex, locPxParamIndex) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.X();
		dF_dXi(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Y();
		dF_dXi(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locDeltaXCrossH_DecayingSource.Z();

		dF_dXi(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier*locQ.X();
		dF_dXi(locFIndex + 1, locPxParamIndex + 1) = locStateSignMultiplier*locQ.Y();
		dF_dXi(locFIndex + 1, locPxParamIndex + 2) = locStateSignMultiplier*locQ.Z();
	}
	else if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle))
	{
		//decaying charged particle in b-field (doesn't include contributions from the particles it is replaced by here)
		if(dDebugLevel > 30)
			cout << "DKinFitter: Calc_dF_Vertex() Section 10; PID = " << locKinFitParticle->Get_PID() << endl;

		TVector3 locQ, locM, locD;
		double locJ;
		Calc_Vertex_Params(locKinFitParticle_DecayingSource, locJ, locQ, locM, locD);

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locBField_DecayingSource = dKinFitUtils->Get_BField(locPosition_DecayingSource);
		TVector3 locH_DecayingSource = locBField_DecayingSource.Unit();
		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
		TVector3 locH = locBField.Unit();
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

		TVector3 locDeltaXCrossH_DecayingSource = locDeltaX_DecayingSource.Cross(locH_DecayingSource);
		TVector3 locDeltaXCrossH_DecayingSource_CrossH = locDeltaXCrossH_DecayingSource.Cross(locH);
		TVector3 locQCrossH = locQ.Cross(locH);

		dF_dXi(locFIndex, locVxParamIndex) -= locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.X();
		dF_dXi(locFIndex, locVxParamIndex + 1) -= locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Y();
		dF_dXi(locFIndex, locVxParamIndex + 2) -= locStateSignMultiplier*locA*locDeltaXCrossH_DecayingSource_CrossH.Z();

		dF_dXi(locFIndex + 1, locVxParamIndex) -= locStateSignMultiplier*locA*locQCrossH.X();
		dF_dXi(locFIndex + 1, locVxParamIndex + 1) -= locStateSignMultiplier*locA*locQCrossH.Y();
		dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locStateSignMultiplier*locA*locQCrossH.Z();

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dXi(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dXi(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dXi(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dXi(locFIndex + 1, locVxParamIndex + 2);
	}
}

void DKinFitter::Calc_dF_Vertex_Decaying_NonAccel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier)
{
	//The input particle is used to (perhaps indirectly) define a decaying particle at the vertex being constrained
		//The decaying particle that this particle directly defines is not accelerating: either neutral or not in a magnetic field
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

	if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && ((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle)))
	{
		//detected charged particle in b-field (can be beam particle)
		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		TVector3 locH = locBField.Unit();

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locStateSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 1) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locStateSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 1) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 11c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locStateSignMultiplier*locDeltaX_DecayingSource.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locStateSignMultiplier*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locStateSignMultiplier*locDeltaX_DecayingSource.X();

			dF_dEta(locFIndex, locVxParamIndex) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dEta(locFIndex, locVxParamIndex + 1) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dEta(locFIndex, locVxParamIndex + 2) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = -1.0*locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag() && (locKinFitParticle->Get_Mass() < 0.00001))
	{
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		TVector3 locPiCrossDeltaXX = locMomentum.Cross(locDeltaX_DecayingSource);
		TVector3 locDeltaXCrossDeltaXX = locDeltaX.Cross(locDeltaX_DecayingSource);
		double locDeltaXiMagSq = locDeltaX.Mag2();

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.X()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.Y()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq - locMomentum.Y()*locDeltaX_DecayingSource.Z()/locDeltaX.Y();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq + locMomentum.Z()*locDeltaX_DecayingSource.Y()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq + locMomentum.X()*locDeltaX_DecayingSource.Z()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq;
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq - locMomentum.Z()*locDeltaX_DecayingSource.X()/locDeltaX.Z();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.X()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.Z()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq - locMomentum.Y()*locDeltaX_DecayingSource.Z()/locDeltaX.Y();
			dF_dEta(locFIndex, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.X()/locDeltaXiMagSq + locMomentum.Z()*locDeltaX_DecayingSource.Y()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq - locMomentum.X()*locDeltaX_DecayingSource.Y()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq + locMomentum.Y()*locDeltaX_DecayingSource.X()/locDeltaX.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq;
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 12c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.Y()/locMomentum.Mag2();
			dF_dEta(locFIndex + 1, locEParamIndex) = locStateSignMultiplier*locEnergy*locPiCrossDeltaXX.Z()/locMomentum.Mag2();

			dF_dEta(locFIndex, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq + locMomentum.X()*locDeltaX_DecayingSource.Z()/locDeltaX.X();
			dF_dEta(locFIndex, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq;
			dF_dEta(locFIndex, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Y()/locDeltaXiMagSq - locMomentum.Z()*locDeltaX_DecayingSource.X()/locDeltaX.Z();
			dF_dEta(locFIndex + 1, locVxParamIndex) = locStateSignMultiplier*locMomentum.X()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq - locMomentum.X()*locDeltaX_DecayingSource.Y()/locDeltaX.X();
			dF_dEta(locFIndex + 1, locVxParamIndex + 1) = locStateSignMultiplier*locMomentum.Y()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq + locMomentum.Y()*locDeltaX_DecayingSource.X()/locDeltaX.Y();
			dF_dEta(locFIndex + 1, locVxParamIndex + 2) = locStateSignMultiplier*locMomentum.Z()*locDeltaXCrossDeltaXX.Z()/locDeltaXiMagSq;
		}

		dF_dXi(locFIndex, locCommonVxParamIndex) -= dF_dEta(locFIndex, locVxParamIndex);
		dF_dXi(locFIndex, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex, locVxParamIndex + 1);
		dF_dXi(locFIndex, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex, locVxParamIndex + 2);

		dF_dXi(locFIndex + 1, locCommonVxParamIndex) -= dF_dEta(locFIndex + 1, locVxParamIndex);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 1) -= dF_dEta(locFIndex + 1, locVxParamIndex + 1);
		dF_dXi(locFIndex + 1, locCommonVxParamIndex + 2) -= dF_dEta(locFIndex + 1, locVxParamIndex + 2);
	}
	else if(locKinFitParticle->Get_IsNeutralShowerFlag() && (locKinFitParticle->Get_Mass() >= 0.00001))
	{
		return;
	}
	else if((locKinFitParticleType == d_DetectedParticle) || (locKinFitParticleType == d_BeamParticle))
	{
		//detected/beam, non-accelerating particle
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex + 1, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 13c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dEta(locFIndex, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dEta(locFIndex, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
			dF_dEta(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dEta(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
	}
	else if((locKinFitParticleType == d_MissingParticle) || ((locKinFitParticleType == d_DecayingParticle) && (locKinFitParticle->Get_PxParamIndex() >= 0)))
	{
		//missing, or open-ended decaying particle
		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex + 1, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locPxParamIndex + 1) = locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = -1.0*locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 14c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locPxParamIndex) = -1.0*locDeltaX_DecayingSource.Z();
			dF_dXi(locFIndex, locPxParamIndex + 2) = locDeltaX_DecayingSource.X();
			dF_dXi(locFIndex + 1, locPxParamIndex) = locDeltaX_DecayingSource.Y();
			dF_dXi(locFIndex + 1, locPxParamIndex + 1) = -1.0*locDeltaX_DecayingSource.X();
		}
	}
	else if((locCharge != 0) && dKinFitUtils->Get_IsBFieldNearBeamline() && (locKinFitParticleType == d_DecayingParticle))
	{
		//decaying charged particle in b-field (doesn't include contributions from the particles it is replaced by here)
		TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
		double locA = -0.00299792458*(double(locCharge))*locBField.Mag();
		TVector3 locH = locBField.Unit();

		TVector3 locPosition_DecayingSource = locKinFitParticle_DecayingSource->Get_Position();
		TVector3 locCommonVertex_DecayingSource = locKinFitParticle_DecayingSource->Get_CommonVertex();
		TVector3 locDeltaX_DecayingSource = locCommonVertex_DecayingSource - locPosition_DecayingSource;

		unsigned short int locVertexConstraintFlag = locKinFitParticle->Get_VertexConstraintFlag();
		if(locVertexConstraintFlag == 1) //1 & 2 //pz is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15a; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 1) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
		}
		else if(locVertexConstraintFlag == 2) //1 & 3 //py is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15b; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.Y()*locH.Y() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 1) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.X();
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.X();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
		}
		else //2 & 3 //px is largest
		{
			if(dDebugLevel > 30)
				cout << "DKinFitter: Calc_dF_Vertex() Section 15c; PID = " << locKinFitParticle->Get_PID() << endl;

			dF_dXi(locFIndex, locVxParamIndex) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Y();
			dF_dXi(locFIndex, locVxParamIndex + 1) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Z()*locH.Z());
			dF_dXi(locFIndex, locVxParamIndex + 2) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Z()*locH.Y();
			dF_dXi(locFIndex + 1, locVxParamIndex) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.X()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 1) += locA*locStateSignMultiplier*locDeltaX_DecayingSource.Y()*locH.Z();
			dF_dXi(locFIndex + 1, locVxParamIndex + 2) -= locA*locStateSignMultiplier*(locDeltaX_DecayingSource.X()*locH.X() + locDeltaX_DecayingSource.Y()*locH.Y());
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

	TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
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

TVector3 DKinFitter::Calc_VertexParams_P4DerivedAtCommonVertex(const DKinFitParticle* locKinFitParticle)
{
	int locCharge = locKinFitParticle->Get_Charge();
	TVector3 locPosition = locKinFitParticle->Get_Position(); //x_0X
	TVector3 locCommonVertex = locKinFitParticle->Get_CommonVertex(); //x_X
	TVector3 locDeltaX = locCommonVertex - locPosition;
	TVector3 locP0X = locKinFitParticle->Get_Momentum(); //at x_0X

	TVector3 locBField = dKinFitUtils->Get_BField(locPosition);
	TVector3 locH = locBField.Unit();
	double locA = -0.00299792458*(double(locCharge))*locBField.Mag();

	TVector3 locDeltaXCrossH = locDeltaX.Cross(locH);
	TVector3 locPX = locP0X - locA*locDeltaXCrossH; //at common vertex x_X

	TVector3 locP0XCrossH = locP0X.Cross(locH);
	double locPXDotH = locPX.Dot(locH);
	double locP0XCrossHMagSq = locP0XCrossH.Mag2();

	double locDeltaXDotH = locDeltaX.Dot(locH);
	double locPXDotDeltaX = locPX.Dot(locDeltaX);

	double locK = locPXDotDeltaX - locPXDotH*locDeltaXDotH;
	double locJ = locA*locK/locP0XCrossHMagSq;

	double locC = locPXDotH/(locP0XCrossHMagSq*sqrt(1.0 - locJ*locJ));
	TVector3 locD = locC*(locPX - locPXDotH*locH) - locH;

	TVector3 locR = 2.0*locJ*locC*locP0XCrossH + locD;
	return locR;
}

/************************************************************************ UPDATE ***********************************************************************/

void DKinFitter::Update_ParticleParams(void)
{
	// update common vertices
	set<DKinFitConstraint_Vertex*> locVertexConstraints = Get_Constraints<DKinFitConstraint_Vertex>();
	set<DKinFitConstraint_Vertex*>::iterator locVertexIterator = locVertexConstraints.begin();
	for(; locVertexIterator != locVertexConstraints.end(); ++locVertexIterator)
	{
		DKinFitConstraint_Vertex* locConstraint = *locVertexIterator;
		int locParamIndex = locConstraint->Get_CommonVxParamIndex();
		TVector3 locVertex(dXi(locParamIndex, 0), dXi(locParamIndex + 1, 0), dXi(locParamIndex + 2, 0));
		locConstraint->Set_CommonVertex(locVertex);
	}

	// update common times
	set<DKinFitConstraint_Spacetime*> locSpacetimeConstraints = Get_Constraints<DKinFitConstraint_Spacetime>();
	set<DKinFitConstraint_Spacetime*>::iterator locSpacetimeIterator = locSpacetimeConstraints.begin();
	for(; locSpacetimeIterator != locSpacetimeConstraints.end(); ++locSpacetimeIterator)
	{
		DKinFitConstraint_Spacetime* locConstraint = *locSpacetimeIterator;
		int locParamIndex = locConstraint->Get_CommonTParamIndex();
		locConstraint->Set_CommonTime(dXi(locParamIndex, 0));
	}

	// update particle information
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if(locKinFitParticleType == d_TargetParticle)
			continue; //nothing to update

		bool locIsUnknownParticleFlag = ((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle));
		TMatrixD& locSourceMatrix = locIsUnknownParticleFlag ? dXi : dEta;

		TVector3 locMomentum = locKinFitParticle->Get_Momentum();
		TVector3 locPreviousMomentum = locMomentum;
		double locPreviousPathLength = locKinFitParticle->Get_PathLength();

		int locParamIndex = locKinFitParticle->Get_PxParamIndex();
		if(locParamIndex >= 0)
		{
			locMomentum = TVector3(locSourceMatrix(locParamIndex, 0), locSourceMatrix(locParamIndex + 1, 0), locSourceMatrix(locParamIndex + 2, 0));
			locKinFitParticle->Set_Momentum(locMomentum);
		}

		//vertex
		locParamIndex = locKinFitParticle->Get_VxParamIndex();
		TVector3 locDeltaX;
		if(locParamIndex >= 0)
		{
			TVector3 locPosition(locSourceMatrix(locParamIndex, 0), locSourceMatrix(locParamIndex + 1, 0), locSourceMatrix(locParamIndex + 2, 0));
			locDeltaX = locPosition - locKinFitParticle->Get_Position();
			locKinFitParticle->Set_Position(locPosition);
			if(locKinFitParticle->Get_CommonVxParamIndex() < 0)
				locKinFitParticle->Set_CommonVertex(locPosition);
		}

		//energy
		locParamIndex = locKinFitParticle->Get_EParamIndex();
		if(locParamIndex >= 0) //neutral shower: set momentum also //must be after Vx & common vertex are set
		{
			double locE = locSourceMatrix(locParamIndex, 0);
			locKinFitParticle->Set_ShowerEnergy(locE);
			double locPMag = sqrt(locE*locE - locKinFitParticle->Get_Mass()*locKinFitParticle->Get_Mass());
			locMomentum = locKinFitParticle->Get_Position() - locKinFitParticle->Get_CommonVertex();
			locMomentum.SetMag(locPMag);
			locKinFitParticle->Set_Momentum(locMomentum);
		}

		//time
		locParamIndex = locKinFitParticle->Get_TParamIndex();
		if(locParamIndex >= 0)
		{
			double locTime = locSourceMatrix(locParamIndex, 0);
			locKinFitParticle->Set_Time(locTime);
			if(locKinFitParticle->Get_CommonTParamIndex() < 0)
				locKinFitParticle->Set_CommonTime(locTime);
		}
		else if(locKinFitParticle->Get_VxParamIndex() >= 0)
		{
			//vertex changed, but time is not a fit parameter, so must manually change time
			double locTime = locKinFitParticle->Get_Time();
			double locNewPathLength = locPreviousPathLength;

			if((locKinFitParticle->Get_Charge() != 0) && dKinFitUtils->Get_IsBFieldNearBeamline()) //in b-field & charged
			{
				TVector3 locH = dKinFitUtils->Get_BField(locKinFitParticle->Get_Position()).Unit();
				double locDeltaXDotH = locDeltaX.Dot(locH);
				double locPDotH = locMomentum.Dot(locH);
				locTime += locDeltaXDotH*locP4.E()/(29.9792458*locPDotH);
			}
			else //non-accelerating
			{
				double locP3Angle = locPreviousMomentum.Angle(locMomentum);
				double locDistanceAtSamePathLength = locPreviousPathLength*sqrt(2.0 - 2.0*cos(locP4Angle));
				double locAlpha = 0.5*(TMath::Pi() - locP4Angle);
				double locF = TMath::Pi() - locAlpha;
				double locDeltaXP3Angle = locDeltaX.Angle(locMomentum);
				double locDeltaXDAngle = TMath::Pi() - locF - locDeltaXP3Angle;
				double locDeltaPathLength = sqrt(locD*locD + locDeltaX.Mag2() - 2*locDeltaX.Mag()*locD*(1.0 - cos(locDeltaXDAngle))); //l3
				double locDeltaT = locDeltaPathLength*locP4.E()/(29.9792458*locMomentum.Mag());

				bool locFurtherAlongTrackFlag = (locDeltaX.Dot(locMomentum) > 0.0);
				locTime += locFurtherAlongTrackFlag ? locDeltaT : -1.0*locDeltaT;
				locNewPathLength += locFurtherAlongTrackFlag ? -1.0*locDeltaPathLength : locDeltaPathLength;
				//double locDeltaXDotP = locDeltaX.Dot(locP4.Vect());
				//locTime += locDeltaXDotP*locP4.E()/(29.9792458*locP4.Vect().Mag2());
			}

			locKinFitParticle->Set_Time(locTime);
			locKinFitParticle->Set_PathLength(locNewPathLength);
		}
	}

	// calc non-unknown decaying particle momentum (derived from other particles)
		//must do last because assumes all other particle p3's are updated
	for(locParticleIterator = dKinFitParticles.begin(); locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue;
		locKinFitParticle->Set_Momentum(dKinFitUtils->Calc_DecayingP4_ByPosition(locKinFitParticle, true).Vect());
	}
}

/************************************************************************ FINAL ************************************************************************/

void DKinFitter::Calc_Pulls(void)
{
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_TargetParticle))
			continue;

		map<DKinFitPullType, double> locParticlePulls;
		if(dDebugLevel >= 50)
		{
			cout << "pulls: PID = " << locKinFitParticle->Get_PID() << endl;
			cout << "e, px, vx, t param indices = " << locKinFitParticle->Get_EParamIndex() << ", " << locKinFitParticle->Get_PxParamIndex() << ", " << locKinFitParticle->Get_VxParamIndex() << ", " << locKinFitParticle->Get_TParamIndex() << endl;
		}

		int locParamIndex = locKinFitParticle->Get_EParamIndex();
		if(locParamIndex >= 0) //E
		{
			double locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_EPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		locParamIndex = locKinFitParticle->Get_PxParamIndex();
		if(locParamIndex >= 0) //px, py, pz
		{
			double locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
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
			double locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
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
			double locDenominator = sqrt(fabs(dVY(locParamIndex, locParamIndex) - (*dVEta)(locParamIndex, locParamIndex)));
			locParticlePulls[d_TPull] = (locDenominator > 0.0) ? dEpsilon(locParamIndex, 0)/locDenominator : std::numeric_limits<double>::quiet_NaN();
		}

		if(!locParticlePulls.empty())
			dPulls[locKinFitParticle] = locParticlePulls;
	}

	if(dDebugLevel > 20)
	{
		cout << "DKinFitter: dEpsilon: " << endl;
		Print_Matrix(dEpsilon);
		cout << "DKinFitter: dVY: " << endl;
		Print_Matrix(dVY);
		cout << "DKinFitter: Pulls: " << endl;
		map<DKinFitParticle*, map<DKinFitPullType, double> >::iterator locIterator;
		for(locIterator = dPulls.begin(); locIterator != dPulls.end(); ++locIterator)
		{
			map<DKinFitPullType, double>& locParticlePulls = locIterator->second;
			DKinFitParticle* locKinFitParticle = locIterator->first;
			TVector3 locMomentum = locKinFitParticle->Get_Momentum();
			cout << "particle PID, p3 = " << locKinFitParticle->Get_PID() << ", " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ":" << endl;
			for(size_t loc_i = 0; loc_i < 8; ++loc_i)
			{
				if(locParticlePulls.find((DKinFitPullType)loc_i) != locParticlePulls.end())
					cout << locParticlePulls[(DKinFitPullType)loc_i] << ", ";
			}
			cout << endl;
		}
	}
}

void DKinFitter::Set_FinalTrackInfo(void)
{
	// first update the covariance matrices of each particle with the fit results (prior to any propagation)
		//correlations between fit and unfit measured parameters: 
			//assume that the correlations remain unchanged: the changes in the parameters and their uncertainties should be small wrst their values
	set<DKinFitParticle*>::iterator locParticleIterator = dKinFitParticles.begin();
	for(; locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if(locKinFitParticleType == d_TargetParticle)
			continue;

		//Get covariance matrix
		bool locReconstructedParticleFlag = ((locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle));
		TMatrixDSym& locKinFitMatrix = locReconstructedParticleFlag ? *dVXi : *dVEta;
		if(locReconstructedParticleFlag) //Brand new particle: Set the covariance matrix from scratch
		{
			//Particle had none: Make a new one
			TMatrixDSym* locCovarianceMatrix = dKinFitUtils->Get_MatrixDSymResource();
			locCovarianceMatrix->ResizeTo(7, 7);
			locCovarianceMatrix->Zero();
			locKinFitParticle->Set_CovarianceMatrix(locCovarianceMatrix);
		}
		//need to update the values, so don't call Get_CovarianceMatrix() function (returns const matrix) //saves memory this way
		TMatrixDSym& locCovarianceMatrix = *(locKinFitParticle->dCovarianceMatrix);

		int locPxParamIndex = locKinFitParticle->Get_PxParamIndex();
		int locVxParamIndex = locKinFitParticle->Get_VxParamIndex();
		int locTParamIndex = locKinFitParticle->Get_TParamIndex();
		int locEParamIndex = locKinFitParticle->Get_EParamIndex();

		int locCovMatrixEParamIndex = locKinFitParticle->Get_CovMatrixEParamIndex();
		int locCovMatrixPxParamIndex = locKinFitParticle->Get_CovMatrixPxParamIndex();
		int locCovMatrixVxParamIndex = locKinFitParticle->Get_CovMatrixVxParamIndex();
		int locCovMatrixTParamIndex = locKinFitParticle->Get_CovMatrixTParamIndex();

		if(dDebugLevel >= 50)
		{
			cout << "SETTING FINAL TRACK INFO: PID = " << locKinFitParticle->Get_PID() << endl;
			cout << "E, px, vx, t param indices = " << locEParamIndex << ", " << locPxParamIndex << ", " << locVxParamIndex << ", " << locTParamIndex << endl;
			cout << "E, px, vx, t cov matrix indices = " << locCovMatrixEParamIndex << ", " << locCovMatrixPxParamIndex << ", " << locCovMatrixVxParamIndex << ", " << locCovMatrixTParamIndex << endl;
			cout << "sizes = " << locCovarianceMatrix.GetNrows() << ", " << locCovarianceMatrix.GetNcols() << ", " << locKinFitMatrix.GetNrows() << ", " << locKinFitMatrix.GetNcols() << endl;
		}

		//decaying particles with momentum derived from other particles
		if((locKinFitParticleType == d_DecayingParticle) && (locPxParamIndex < 0))
		{
			//decaying particle: the momentum is derived from the other particles
			//if the decaying particle is not a constraining particle (e.g. vertex no-constrain):
			//look for particles that define the momentum that were NOT in the fit. They're uncertainties are needed, but are not in the input matrices

			set<DKinFitParticle*> locDerivedFromParticles = locKinFitParticle->Get_FromAllParticles();
			set<DKinFitParticle*>::iterator locFromIterator = locDerivedFromParticles.begin();
			map<DKinFitParticle*, int> locAdditionalPxParamIndices;
			for(; locFromIterator != locDerivedFromParticles.end(); ++locFromIterator)
			{
				DKinFitParticleType locKinFitParticleType = (*locFromIterator)->Get_KinFitParticleType();
				if((locKinFitParticleType == d_DecayingParticle) || (locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle))
					continue;
				if((*locFromIterator)->Get_PxParamIndex() >= 0)
					continue; //uncertainties are already in matrices: ignore
				int locNewPxParamIndex = dNumEta + dNumXi + 3*locAdditionalPxParamIndices.size();
				locAdditionalPxParamIndices[*locFromIterator] = locNewPxParamIndex;
			}

			TMatrixD locJacobian(7, dNumEta + dNumXi + 3*locAdditionalPxParamIndices.size());
			locJacobian.Zero();

			bool locP3DerivedAtProductionVertexFlag = !locKinFitParticle->Get_FromInitialState().empty(); //else decay vertex
			bool locP3DerivedAtPositionFlag = (locP3DerivedAtProductionVertexFlag == locKinFitParticle->Get_VertexP4AtProductionVertex());
			dKinFitUtils->Calc_DecayingParticleJacobian(locKinFitParticle, locP3DerivedAtPositionFlag, 1.0, dNumEta, locAdditionalPxParamIndices, locJacobian);

			//set vertex & time terms
			if(locVxParamIndex >= 0)
			{
				locJacobian(3, locVxParamIndex + dNumEta) = 1.0;
				locJacobian(4, locVxParamIndex + dNumEta + 1) = 1.0;
				locJacobian(5, locVxParamIndex + dNumEta + 2) = 1.0;
			}
			int locTParamIndex = locKinFitParticle->Get_TParamIndex();
			if(locTParamIndex >= 0)
				locJacobian(6, locTParamIndex + dNumEta) = 1.0;

			if(dDebugLevel >= 50)
			{
				cout << "JACOBIAN MATRIX (enclosed decaying particle):" << endl;
				Print_Matrix(locJacobian);
			}

			//build matrix, transform errors
			if(locAdditionalPxParamIndices.empty())
			{
				TMatrixDSym locTempCovarianceMatrix = *dV;
				locCovarianceMatrix = locTempCovarianceMatrix.Similarity(locJacobian);
			}
			else //need to expand error matrix
			{
				TMatrixDSym locTempCovarianceMatrix(dNumEta + dNumXi + 3*locAdditionalPxParamIndices.size());
				locTempCovarianceMatrix.SetSub(0, *dV); //insert dV

				//insert p3 covariance matrices of additional particles
				map<DKinFitParticle*, int>::iterator locPxParamIterator = locAdditionalPxParamIndices.begin();
				for(; locPxParamIterator != locAdditionalPxParamIndices.end(); ++locPxParamIterator)
				{
					int locPxParamIndex = locPxParamIterator->second;
					TMatrixDSym locAdditionalCovMatrix(3);
					locAdditionalCovMatrix = locPxParamIterator->first->Get_CovarianceMatrix()->GetSub(0, 2, locAdditionalCovMatrix);
					locTempCovarianceMatrix.SetSub(locPxParamIndex, locAdditionalCovMatrix);
				}

				if(dDebugLevel >= 50)
				{
					cout << "FULL ERROR MATRIX (for calculated cov for enclosed decaying particle):" << endl;
					Print_Matrix(locTempCovarianceMatrix);
				}

				//transform errors
				locCovarianceMatrix = locTempCovarianceMatrix.Similarity(locJacobian);
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
			double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/locDenominator : 0.0;
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
				double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
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
			double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locEParamIndex, locEParamIndex))/locDenominator : 0.0;
			locCovarianceMatrix(locCovMatrixEParamIndex, locCovMatrixTParamIndex) *= locUncertaintyRatio;
			locCovarianceMatrix(locCovMatrixTParamIndex, locCovMatrixEParamIndex) *= locUncertaintyRatio;
		}
		else if(!locReconstructedParticleFlag && (locEParamIndex < 0) && (locTParamIndex >= 0) && (locCovMatrixEParamIndex >= 0)) //only T included in the fit //E may not be in the covariance matrix!!
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
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
				double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/locDenominator : 0.0;
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
				double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
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
				double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locPxParamIndex + loc_j, locPxParamIndex + loc_j))/locDenominator : 0.0;
				locCovarianceMatrix(locCovMatrixPxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixPxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locPxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
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
				double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locVxParamIndex + loc_j, locVxParamIndex + loc_j))/locDenominator : 0.0;
				locCovarianceMatrix(locCovMatrixVxParamIndex + loc_j, locCovMatrixTParamIndex + 0) *= locUncertaintyRatio;
				locCovarianceMatrix(locCovMatrixTParamIndex + 0, locCovMatrixVxParamIndex + loc_j) *= locUncertaintyRatio;
			}
		}
		else if(!locReconstructedParticleFlag && (locVxParamIndex < 0) && (locTParamIndex >= 0)) //only T included in the fit
		{
			double locDenominator = sqrt(dVY(locTParamIndex, locTParamIndex));
			double locUncertaintyRatio = (fabs(locDenominator) > 0.0) ? sqrt(locKinFitMatrix(locTParamIndex, locTParamIndex))/locDenominator : 0.0;
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
	for(locParticleIterator = dKinFitParticles.begin(); locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();
		if((locKinFitParticleType == d_TargetParticle) || (locKinFitParticleType == d_MissingParticle) || (locKinFitParticleType == d_DecayingParticle))
			continue; // particle properties already defined at the fit vertex

		if(!locKinFitParticle->Get_FitCommonVertexFlag())
			continue; // no distance over which to propagate

		pair<double, double> locPathLengthPair;
		TMatrixDSym& locCovarianceMatrix = *(locKinFitParticle->dCovarianceMatrix);

		TVector3 locMomentum;
		TLorentzVector locSpacetimeVertex;
		if(!dKinFitUtils->Propagate_TrackInfoToCommonVertex(locKinFitParticle, dVXi, locMomentum, locSpacetimeVertex, locPathLengthPair, locCovarianceMatrix))
			continue; // info not propagated

		if(dDebugLevel >= 50)
		{
			cout << "PROPAGATED FINAL TRACK INFO: PID = " << locKinFitParticle->Get_PID() << endl;
			cout << "p_xyz, v_xyzt = " << locMomentum.Px() << ", " << locMomentum.Py() << ", " << locMomentum.Pz() << ", " << locSpacetimeVertex.X() << ", " << locSpacetimeVertex.Y() << ", " << locSpacetimeVertex.Z() << ", " << locSpacetimeVertex.T() << endl;
			cout << "common v_xyzt = " << locKinFitParticle->Get_CommonVertex().X() << ", " << locKinFitParticle->Get_CommonVertex().Y() << ", " << locKinFitParticle->Get_CommonVertex().Z() << ", " << locKinFitParticle->Get_CommonTime() << endl;
			cout << "path length & uncert = " << locPathLengthPair.first << ", " << locPathLengthPair.second << endl;
			cout << "sizes = " << locCovarianceMatrix.GetNrows() << ", " << locCovarianceMatrix.GetNcols() << endl;
		}

		//no need to set the covariance matrix: already updated (passed in a reference to it)
		locKinFitParticle->Set_Momentum(locMomentum);
		if(locKinFitParticle->Get_IsNeutralShowerFlag())
			locKinFitParticle->Set_CommonSpacetimeVertex(locSpacetimeVertex);
		else
			locKinFitParticle->Set_SpacetimeVertex(locSpacetimeVertex);
		locKinFitParticle->Set_PathLength(locPathLengthPair.first);
		locKinFitParticle->Set_PathLengthUncertainty(locPathLengthPair.second);
	}

	//calculate the path length of decaying particles involved in 2 vertex fits
	for(locParticleIterator = dKinFitParticles.begin(); locParticleIterator != dKinFitParticles.end(); ++locParticleIterator)
	{
		DKinFitParticle* locKinFitParticle = *locParticleIterator;
		DKinFitParticleType locKinFitParticleType = locKinFitParticle->Get_KinFitParticleType();

		if((locKinFitParticleType != d_DecayingParticle) || !locKinFitParticle->Get_FitCommonVertexFlag())
			continue;

		pair<double, double> locPathLengthPair;
		const TMatrixDSym& locCovarianceMatrix = *(locKinFitParticle->Get_CovarianceMatrix());
		if(dKinFitUtils->Calc_PathLength(locKinFitParticle, dVXi, locCovarianceMatrix, locPathLengthPair))
		{
			locKinFitParticle->Set_PathLength(locPathLengthPair.first);
			locKinFitParticle->Set_PathLengthUncertainty(locPathLengthPair.second);
		}
	}
}

