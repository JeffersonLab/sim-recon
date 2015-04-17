#ifndef _DKinFitter_
#define _DKinFitter_

#include <deque>
#include <utility>
#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <limits>

#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TDecompLU.h"

#include "DKinFitConstraints.h"
#include "DKinFitParticle.h"

using namespace std;

enum DKinFitStatus
{
	d_KinFitSuccessful = 0,
	d_KinFitFailedSetup,
	d_KinFitFailedInversion,
	d_KinFitTooManyIterations
};

class DKinFitter //purely virtual: cannot directly instantiate class, can only inherit from it
{
	protected:

		DKinFitter(void);
		virtual ~DKinFitter(void);

	public:

		virtual void Reset_NewEvent(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!
		virtual void Reset_NewFit(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!

		void Preallocate_MatrixMemory(void);

		//CREATE PARTICLES
			//If multiple constraints, it is EXTREMELY CRITICAL that only one DKinFitParticle be created per particle, so that the particles are correctly linked across constraints!!
			//data is not registered with the fitter (memory is though)
		const DKinFitParticle* Make_DecayingParticle(int locPID, int locCharge, double locMass);
		const DKinFitParticle* Make_MissingParticle(int locPID, int locCharge, double locMass);
		const DKinFitParticle* Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		const DKinFitParticle* Make_TargetParticle(int locPID, int locCharge, double locMass);
		const DKinFitParticle* Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		const DKinFitParticle* Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix);

		//INPUT RF TIME
		void Set_RFTime(double locRFTime, double locRFUncertainty, const DKinFitParticle* locRFMatchedBeamParticle);

		//CREATE CONSTRAINTS
			//data is not registered with the fitter (memory is though)
		DKinFitConstraint_Vertex* Make_VertexConstraint(const deque<const DKinFitParticle*>& locFinalParticles, TVector3 locVertexGuess);
		DKinFitConstraint_Vertex* Make_VertexConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, TVector3 locVertexGuess);
		DKinFitConstraint_Spacetime* Make_SpacetimeConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locUseRFTimeFlag, TVector3 locVertexGuess, double locCommonTimeGuess);
		//note, below locConstrainInitialParticleMassFlag is ignored if the parent particle is beam/detected/open-ended-decaying (enforce p4 constraint instead)
		DKinFitConstraint_P4* Make_P4Constraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locConstrainInitialParticleMassFlag = true);

		//CLONE CONSTRAINTS
			//if need to modify a constraint without disrupting the original: note that particles aren't cloned!
		DKinFitConstraint_P4* Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint);
		DKinFitConstraint_Vertex* Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint);
		DKinFitConstraint_Spacetime* Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint);

		//SORT CONSTRAINTS
			//used to determine which order to perform partial fits to get vertex guesses prior to the full fit
		bool Sort_Constraints(const deque<DKinFitConstraint*>& locOriginalConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints);

		//SET CONSTRAINTS FOR FIT
		void Set_Constraint(const DKinFitConstraint_P4* locConstraint);
		void Set_Constraint(const DKinFitConstraint_Vertex* locConstraint);
		void Set_Constraint(const DKinFitConstraint_Spacetime* locConstraint);

		//FIT!
		virtual bool Fit_Reaction(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!

		//GET STATUS (Fit_Reaction returns true/false for success/fail, but this provides more details on the failures
		DKinFitStatus Get_KinFitStatus(void) const{return dKinFitStatus;}

		//GET/SET CONTROL VARIABLES
		void Set_ConvergenceChiSqDiff(double locConvergenceChiSqDiff){dConvergenceChiSqDiff = locConvergenceChiSqDiff;}
		double Get_ConvergenceChiSqDiff(void) const{return dConvergenceChiSqDiff;}
		void Set_ConvergenceChiSqDiff_LastResort(double locConvergenceChiSqDiff){dConvergenceChiSqDiff_LastResort = locConvergenceChiSqDiff;}
		double Get_ConvergenceChiSqDiff_LastResort(void) const{return dConvergenceChiSqDiff_LastResort;}
		void Set_LinkVerticesFlag(bool locLinkVerticesFlag){dLinkVerticesFlag = locLinkVerticesFlag;}
		bool Get_LinkVerticesFlag(void) const{return dLinkVerticesFlag;}
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}
		int Get_DebugLevel(void) const{return dDebugLevel;}
		void Set_MaxNumIterations(unsigned int locMaxNumIterations){dMaxNumIterations = locMaxNumIterations;}
		unsigned int Get_MaxNumIterations(void) const{return dMaxNumIterations;}
		void Set_MaxKinFitParticlePoolSize(size_t locMaxKinFitParticlePoolSize){dMaxKinFitParticlePoolSize = locMaxKinFitParticlePoolSize;}
		size_t Get_MaxKinFitParticlePoolSize(void) const{return dMaxKinFitParticlePoolSize;}
		void Set_MaxKinFitConstraintVertexPoolSize(size_t locMaxKinFitConstraintVertexPoolSize){dMaxKinFitConstraintVertexPoolSize = locMaxKinFitConstraintVertexPoolSize;}
		size_t Get_MaxKinFitConstraintVertexPoolSize(void) const{return dMaxKinFitConstraintVertexPoolSize;}
		void Set_MaxKinFitConstraintSpacetimePoolSize(size_t locMaxKinFitConstraintSpacetimePoolSize){dMaxKinFitConstraintSpacetimePoolSize = locMaxKinFitConstraintSpacetimePoolSize;}
		size_t Get_MaxKinFitConstraintSpacetimePoolSize(void) const{return dMaxKinFitConstraintSpacetimePoolSize;}
		void Set_MaxKinFitConstraintP4PoolSize(size_t locMaxKinFitConstraintP4PoolSize){dMaxKinFitConstraintP4PoolSize = locMaxKinFitConstraintP4PoolSize;}
		size_t Get_MaxKinFitConstraintP4PoolSize(void) const{return dMaxKinFitConstraintP4PoolSize;}
		void Set_MaxMatrixDSymPoolSize(size_t locMaxMatrixDSymPoolSize){dMaxMatrixDSymPoolSize = locMaxMatrixDSymPoolSize;}
		size_t Get_MaxMatrixDSymPoolSize(void) const{return dMaxMatrixDSymPoolSize;}
		void Set_MaxLargeMatrixDSymPoolSize(size_t locMaxLargeMatrixDSymPoolSize){dMaxLargeMatrixDSymPoolSize = locMaxLargeMatrixDSymPoolSize;}
		size_t Get_MaxLargeMatrixDSymPoolSize(void) const{return dMaxLargeMatrixDSymPoolSize;}

		//GET FIT INFORMATION
		unsigned int Get_NumUnknowns(void) const{return dNumXi;}
		unsigned int Get_NumMeasurables(void) const{return dNumEta;}
		unsigned int Get_NumConstraintEquations(void) const{return dNumF;}

		//GET FIT QUALITY RESULTS
		double Get_ChiSq(void) const{return dChiSq;}
		double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}
		unsigned int Get_NDF(void) const{return dNDF;}
		void Get_Pulls(map<const DKinFitParticle*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;} //key is particle (NULL for rf-bunch), 2nd key is param type

		//GET FIT PARAMETERS
		const DKinFitParticle* Get_OutputKinFitParticle(const DKinFitParticle* locInputKinFitParticle) const;
		const DKinFitParticle* Get_InputKinFitParticle(const DKinFitParticle* locOutputKinFitParticle) const;
		set<const DKinFitParticle*> Get_InputKinFitParticleSet(const DKinFitParticle* locOutputKinFitParticle) const;

		const DKinFitConstraint* Get_OutputKinFitConstraint(const DKinFitConstraint* locInputKinFitConstraint) const;
		const DKinFitConstraint* Get_InputKinFitConstraint(const DKinFitConstraint* locOutputKinFitConstraint) const;
		set<const DKinFitConstraint*> Get_InputKinFitConstraintSet(const DKinFitConstraint* locOutputKinFitConstraint) const;

		void Get_KinFitConstraints(deque<const DKinFitConstraint*>& locKinFitConstraints) const;
		void Get_KinFitParticles(deque<const DKinFitParticle*>& locKinFitParticles) const;

		double Get_RFTime(void) const{return dRFTime;}

		//GET UNCERTAINTIES
		const TMatrixDSym* Get_VEta(void) {return dVEta;}
		const TMatrixDSym* Get_VXi(void) {return dVXi;}
		const TMatrixDSym* Get_V(void) {return dV;}

		//GET CURRENT POOL SIZES
		size_t Get_KinFitParticlePoolSize(void) const{return dKinFitParticlePool_All.size();};
		size_t Get_KinFitConstraintVertexPoolSize(void) const{return dKinFitConstraintVertexPool_All.size();};
		size_t Get_KinFitConstraintSpacetimePoolSize(void) const{return dKinFitConstraintSpacetimePool_All.size();};
		size_t Get_KinFitConstraintP4PoolSize(void) const{return dKinFitConstraintP4Pool_All.size();};
		size_t Get_MatrixDSymPoolSize(void) const{return dMatrixDSymPool_All.size();};
		size_t Get_LargeMatrixDSymPoolSize(void) const{return dLargeMatrixDSymPool_All.size();};

		virtual bool Get_IsBFieldNearBeamline(void) const = 0;

		bool Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixDSym& locCovarianceMatrix) const;
		bool Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixDSym& locCovarianceMatrix, pair<double, double>& locPathLengthPair) const;

	protected:
		virtual TVector3 Get_BField(const TVector3& locPosition) const = 0; //must return in units of Tesla!!
		TMatrixDSym* Get_MatrixDSymResource(void);

	private:

		DKinFitParticle* Clone_KinFitParticle(const DKinFitParticle* locKinFitParticle);
		TMatrixDSym* Clone_MatrixDSym(const TMatrixDSym* locMatrix);

		DKinFitParticle* Get_KinFitParticleResource(void);
		DKinFitConstraint_Vertex* Get_KinFitConstraintVertexResource(void);
		DKinFitConstraint_Spacetime* Get_KinFitConstraintSpacetimeResource(void);
		DKinFitConstraint_P4* Get_KinFitConstraintP4Resource(void);
		TMatrixDSym* Get_LargeMatrixDSymResource(void);

		bool Calc_dS(void);
		bool Calc_dU(void);

		void Clone_ConstraintsForFit(void);
		void Register_ParticlesForFit(void);
		void Set_MatrixSizes(void);
		void Resize_Matrices(void);
		void Zero_Matrices(void);

		// PREPARE CONSTRAINTS
		void Clone_Constraints(const deque<DKinFitConstraint*>& locInputConstraints, deque<DKinFitConstraint*>& locClonedConstraints);
		bool Prepare_Constraint(DKinFitConstraint_P4* locConstraint) const;
		bool Prepare_Constraint(DKinFitConstraint_Vertex* locConstraint) const;
		bool Prepare_Constraint(DKinFitConstraint_Spacetime* locConstraint) const;

		// RESOLVE CONSTRAINTS
		bool Resolve_Constraints(void);
		bool Resolve_Constraints(const deque<DKinFitConstraint*>& locConstraints, deque<DKinFitConstraint_VertexBase*>& locSortedVertexConstraints, bool locSortOnlyFlag) const;
		bool Resolve_DecayingParticleSpacetimeLinks(const deque<DKinFitConstraint*>& locKinFitConstraints, deque<DKinFitConstraint_VertexBase*>& locSortedConstraints, bool locSortOnlyFlag) const;
		bool Resolve_P4Constraints(const deque<DKinFitConstraint*>& locKinFitConstraints, bool locSortOnlyFlag) const;
		bool Find_ConstrainableParticles(const deque<DKinFitConstraint_P4*>& locP4Constraints, deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles) const;
		void Constrain_Particle(DKinFitParticle* locParticleToConstrain, DKinFitConstraint_P4* locConstraint, TLorentzVector& locP4) const;
		bool Resolve_P4MassConstraints(const deque<DKinFitConstraint*>& locKinFitConstraints, bool locSortOnlyFlag) const;
		bool Resolve_InclusiveP4(const deque<DKinFitConstraint*>& locKinFitConstraints) const;
		void Mark_AsMissingMassConstraintIfNecessary(DKinFitConstraint_P4* locP4Constraint) const;
		bool Group_Constraints(const deque<DKinFitConstraint_VertexBase*>& locSortedVertexConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints) const;
		void Recycle_Constraints(deque<DKinFitConstraint*>& locConstraints);

		void Fill_InputMatrices(void);
		void Update_ParticleParams(void);
		TLorentzVector Calc_DecayingP4(DKinFitConstraint_P4* locP4Constraint, bool locDecayMomentumFlag) const;

		void Print_Matrix(const TMatrixD& locMatrix) const;
		void Print_ParticleParams(const DKinFitParticle* locKinFitParticle) const;

		void Calc_dF(void);
		void Calc_dF_P4(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, DKinFitConstraint_P4* locKinFitSubConstraint_P4);
		TLorentzVector Calc_dF_MassP4(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, bool locIsInvariantMassConstraint, DKinFitConstraint_P4* locKinFitSubConstraint_P4);
		void Calc_dF_MassDerivs(DKinFitConstraint_P4* locKinFitConstraint_P4, const DKinFitParticle* locKinFitParticle, TLorentzVector locDecayingParticleDerivedP4, bool locInitialStateFlag, bool locIsInvariantMassConstraint, DKinFitConstraint_P4* locKinFitSubConstraint_P4);

		void Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag);
		void Calc_dF_Vertex_NotDecaying(size_t locFIndex, const DKinFitParticle* locKinFitParticle);
		void Calc_dF_Vertex_Decaying_Accel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag);
		void Calc_dF_Vertex_Decaying_NonAccel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, bool locInitialStateFlag, bool locOriginalInitialStateFlag);
		void Calc_Vertex_Params(const DKinFitParticle* locKinFitParticle, double& locJ, TVector3& locQ, TVector3& locM, TVector3& locD);

		void Calc_dF_Time(size_t locFIndex, const DKinFitParticle* locKinFitParticle, bool locUseRFTimeFlag);

		void Calc_Pulls(void);
		void Set_FinalTrackInfo(void);
		void Calc_DecayingParticleJacobian(DKinFitConstraint_P4* locP4Constraint, bool locDecayVertexFlag, TMatrixD& locJacobian) const;

		unsigned int dDebugLevel;

		deque<DKinFitConstraint*> dKinFitConstraints;
		deque<DKinFitParticle*> dKinFitParticles;

		// maps for registering cloning //key = cloned output, key = "input" (which could be a clone!)
			//can't have IO maps: can have multiple clones of a given input object (e.g. multiple fits)
		map<const DKinFitParticle*, const DKinFitParticle*> dKinFitParticleCloningOIMap;
		map<const DKinFitConstraint*, const DKinFitConstraint*> dKinFitConstraintCloningOIMap;

		//RF Time Stuff
		double dRFTime;
		double dRFUncertainty;
		DKinFitParticle* dRFMatchedBeamParticle;
		int dRFTimeParamIndex;

		TMatrixD dXi; //unmeasurable unknowns
		TMatrixD dEta; //observables
		TMatrixD dY; //first approximation of observables (initial measurements)
		TMatrixDSym dVY; //convariance matrix of dY

		TMatrixDSym dS; // the covariance matrix of the uncertainties of whether the constraint equations are satisfied
		TMatrixDSym dS_Inverse;
		TMatrixDSym dU;
		TMatrixDSym dU_Inverse;

		TMatrixD dF; //constraint equations with eta dependence
		TMatrixD dEpsilon; //dY - dEta
		TMatrixD dLambda; //lagrange multipliers of constraint equations
		TMatrixD dLambda_T;

		TMatrixD dF_dEta; //partial derivative of constraint equations wrst the observables
		TMatrixD dF_dEta_T;
		TMatrixD dF_dXi; //partial derivative of constraint equations wrst the unmeasurable unknowns
		TMatrixD dF_dXi_T;

		TMatrixDSym* dVXi; //covariance matrix of dXi
		TMatrixDSym* dVEta; //covariance matrix of dEta
		TMatrixDSym* dV; //full covariance matrix: dVEta at top-left and dVXi at bottom-right (+ the eta, xi covariance)

		unsigned int dNumXi; //num unknowns
		unsigned int dNumEta; //num measurables
		unsigned int dNumF; //num constraint eqs

		double dChiSq;
		unsigned int dNDF;
		double dConfidenceLevel;
		map<const DKinFitParticle*, map<DKinFitPullType, double> > dPulls; //key is particle (NULL for rf-bunch), 2nd key is param type
		//e.g. 2nd dimension can be E, x, y, z, t for neutral showers; px, py, pz, x, y, z, t for charged or neutral tracks; or a subset of these; or t for RF

		deque<DKinFitParticle*> dKinFitParticlePool_All;
		deque<DKinFitParticle*> dKinFitParticlePool_Available;

		deque<DKinFitConstraint_Vertex*> dKinFitConstraintVertexPool_All;
		deque<DKinFitConstraint_Vertex*> dKinFitConstraintVertexPool_Available;

		deque<DKinFitConstraint_Spacetime*> dKinFitConstraintSpacetimePool_All;
		deque<DKinFitConstraint_Spacetime*> dKinFitConstraintSpacetimePool_Available;

		deque<DKinFitConstraint_P4*> dKinFitConstraintP4Pool_All;
		deque<DKinFitConstraint_P4*> dKinFitConstraintP4Pool_Available;

		deque<TMatrixDSym*> dMatrixDSymPool_All;
		deque<TMatrixDSym*> dMatrixDSymPool_Available;

		deque<TMatrixDSym*> dLargeMatrixDSymPool_All;
		deque<TMatrixDSym*> dLargeMatrixDSymPool_Available;

		size_t dMaxKinFitParticlePoolSize;
		size_t dMaxKinFitConstraintVertexPoolSize;
		size_t dMaxKinFitConstraintSpacetimePoolSize;
		size_t dMaxKinFitConstraintP4PoolSize;
		size_t dMaxMatrixDSymPoolSize;
		size_t dMaxLargeMatrixDSymPoolSize;

		unsigned int dMaxNumIterations;
		bool dLinkVerticesFlag;
		DKinFitStatus dKinFitStatus;
		double dConvergenceChiSqDiff;
		double dConvergenceChiSqDiff_LastResort; //if max # iterations hit, use this for final check (sometimes chisq walks (very slightly) forever without any meaningful change in the variables)
};

inline DKinFitConstraint_Vertex* DKinFitter::Make_VertexConstraint(const deque<const DKinFitParticle*>& locFinalParticles, TVector3 locVertexGuess)
{
	deque<const DKinFitParticle*> locInitialParticles;
	return Make_VertexConstraint(locInitialParticles, locFinalParticles, locVertexGuess);
}

inline const DKinFitParticle* DKinFitter::Get_OutputKinFitParticle(const DKinFitParticle* locInputKinFitParticle) const
{
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
	{
		set<const DKinFitParticle*> locInputKinFitParticleSet = Get_InputKinFitParticleSet(dKinFitParticles[loc_i]);
		if(locInputKinFitParticleSet.find(locInputKinFitParticle) != locInputKinFitParticleSet.end())
			return dKinFitParticles[loc_i];
	}
	return NULL;
}

inline const DKinFitParticle* DKinFitter::Get_InputKinFitParticle(const DKinFitParticle* locOutputKinFitParticle) const
{
	const DKinFitParticle* locInputKinFitParticle = locOutputKinFitParticle; //perhaps the argument is the input particle
	const DKinFitParticle* locKinFitParticleToSearchFor = locOutputKinFitParticle;
	while(true)
	{
		map<const DKinFitParticle*, const DKinFitParticle*>::const_iterator locIterator = dKinFitParticleCloningOIMap.find(locKinFitParticleToSearchFor);
		if(locIterator == dKinFitParticleCloningOIMap.end())
			return locInputKinFitParticle; //input particle if none
		locInputKinFitParticle = locIterator->second;
		locKinFitParticleToSearchFor = locInputKinFitParticle;
	}
	return NULL; //impossible to reach here anyway
}

inline set<const DKinFitParticle*> DKinFitter::Get_InputKinFitParticleSet(const DKinFitParticle* locOutputKinFitParticle) const
{
	//set of ALL input particles that led to this output clone (not just the original)
	const DKinFitParticle* locKinFitParticleToSearchFor = locOutputKinFitParticle;
	set<const DKinFitParticle*> locInputKinFitParticleSet;
	locInputKinFitParticleSet.insert(locOutputKinFitParticle); //perhaps the argument is the input particle
	while(true)
	{
		map<const DKinFitParticle*, const DKinFitParticle*>::const_iterator locIterator = dKinFitParticleCloningOIMap.find(locKinFitParticleToSearchFor);
		if(locIterator == dKinFitParticleCloningOIMap.end())
			return locInputKinFitParticleSet; //input particle if none
		locInputKinFitParticleSet.insert(locIterator->second);
		locKinFitParticleToSearchFor = locIterator->second;
	}
	return set<const DKinFitParticle*>(); //impossible to reach here anyway
}

inline const DKinFitConstraint* DKinFitter::Get_OutputKinFitConstraint(const DKinFitConstraint* locInputKinFitConstraint) const
{
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
	{
		set<const DKinFitConstraint*> locInputKinFitConstraintSet = Get_InputKinFitConstraintSet(dKinFitConstraints[loc_i]);
		if(locInputKinFitConstraintSet.find(locInputKinFitConstraint) != locInputKinFitConstraintSet.end())
			return dKinFitConstraints[loc_i];
	}
	return NULL;
}

inline const DKinFitConstraint* DKinFitter::Get_InputKinFitConstraint(const DKinFitConstraint* locOutputKinFitConstraint) const
{
	const DKinFitConstraint* locInputKinFitConstraint = locOutputKinFitConstraint; //perhaps the argument is the input particle
	const DKinFitConstraint* locKinFitConstraintToSearchFor = locOutputKinFitConstraint;
	while(true)
	{
		map<const DKinFitConstraint*, const DKinFitConstraint*>::const_iterator locIterator = dKinFitConstraintCloningOIMap.find(locKinFitConstraintToSearchFor);
		if(locIterator == dKinFitConstraintCloningOIMap.end())
			return locInputKinFitConstraint; //input particle if none
		locInputKinFitConstraint = locIterator->second;
		locKinFitConstraintToSearchFor = locInputKinFitConstraint;
	}
	return NULL; //impossible to reach here anyway
}

inline set<const DKinFitConstraint*> DKinFitter::Get_InputKinFitConstraintSet(const DKinFitConstraint* locOutputKinFitConstraint) const
{
	//set of ALL input particles that led to this output clone (not just the original)
	const DKinFitConstraint* locKinFitConstraintToSearchFor = locOutputKinFitConstraint;
	set<const DKinFitConstraint*> locInputKinFitConstraintSet;
	locInputKinFitConstraintSet.insert(locOutputKinFitConstraint); //perhaps the argument is the input particle
	while(true)
	{
		map<const DKinFitConstraint*, const DKinFitConstraint*>::const_iterator locIterator = dKinFitConstraintCloningOIMap.find(locKinFitConstraintToSearchFor);
		if(locIterator == dKinFitConstraintCloningOIMap.end())
			return locInputKinFitConstraintSet; //input particle if none
		locInputKinFitConstraintSet.insert(locIterator->second);
		locKinFitConstraintToSearchFor = locIterator->second;
	}
	return set<const DKinFitConstraint*>(); //impossible to reach here anyway
}

inline void DKinFitter::Get_KinFitConstraints(deque<const DKinFitConstraint*>& locKinFitConstraints) const
{
	for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
		locKinFitConstraints.push_back(dKinFitConstraints[loc_i]);
}

inline void DKinFitter::Get_KinFitParticles(deque<const DKinFitParticle*>& locKinFitParticles) const
{
	for(size_t loc_i = 0; loc_i < dKinFitParticles.size(); ++loc_i)
		locKinFitParticles.push_back(dKinFitParticles[loc_i]);
}

#endif // _DKinFitter_
