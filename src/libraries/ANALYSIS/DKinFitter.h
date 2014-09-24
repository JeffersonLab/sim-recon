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
		inline DKinFitConstraint_P4* Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint)
		{
			return Clone_KinFitConstraint_P4(locConstraint, false);
		}
		inline DKinFitConstraint_Vertex* Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint)
		{
			return Clone_KinFitConstraint_Vertex(locConstraint, false);
		}
		inline DKinFitConstraint_Spacetime* Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint)
		{
			return Clone_KinFitConstraint_Spacetime(locConstraint, false);
		}

		//SORT CONSTRAINTS
			//used to determine which order to perform partial fits to get vertex guesses prior to the full fit
		bool Sort_Constraints(const deque<DKinFitConstraint*>& locOriginalConstraints, deque<pair<DKinFitConstraint_VertexBase*, set<DKinFitConstraint_P4*> > >& locSortedConstraints);

		//SET CONSTRAINTS FOR FIT
			//registers particles & constraints with the fitter
			//clones input particles & constraints so that the originals aren't modified during fitting. clones are returned
		const DKinFitConstraint_P4* Set_Constraint(const DKinFitConstraint_P4* locConstraint);
		const DKinFitConstraint_Vertex* Set_Constraint(const DKinFitConstraint_Vertex* locConstraint);
		const DKinFitConstraint_Spacetime* Set_Constraint(const DKinFitConstraint_Spacetime* locConstraint);

		//FIT!
		virtual bool Fit_Reaction(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!

		//GET STATUS (Fit_Reaction returns true/false for success/fail, but this provides more details on the failures
		inline DKinFitStatus Get_KinFitStatus(void) const{return dKinFitStatus;}

		//GET/SET CONTROL VARIABLES
		inline void Set_ConvergenceChiSqDiff(double locConvergenceChiSqDiff){dConvergenceChiSqDiff = locConvergenceChiSqDiff;}
		inline double Get_ConvergenceChiSqDiff(void) const{return dConvergenceChiSqDiff;}
		inline void Set_ConvergenceChiSqDiff_LastResort(double locConvergenceChiSqDiff){dConvergenceChiSqDiff_LastResort = locConvergenceChiSqDiff;}
		inline double Get_ConvergenceChiSqDiff_LastResort(void) const{return dConvergenceChiSqDiff_LastResort;}
		inline void Set_LinkVerticesFlag(bool locLinkVerticesFlag){dLinkVerticesFlag = locLinkVerticesFlag;}
		inline bool Get_LinkVerticesFlag(void) const{return dLinkVerticesFlag;}
		inline void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}
		inline int Get_DebugLevel(void) const{return dDebugLevel;}
		inline void Set_MaxNumIterations(unsigned int locMaxNumIterations){dMaxNumIterations = locMaxNumIterations;}
		inline unsigned int Get_MaxNumIterations(void) const{return dMaxNumIterations;}
		inline void Set_MaxKinFitParticlePoolSize(size_t locMaxKinFitParticlePoolSize){dMaxKinFitParticlePoolSize = locMaxKinFitParticlePoolSize;}
		inline size_t Get_MaxKinFitParticlePoolSize(void) const{return dMaxKinFitParticlePoolSize;}
		inline void Set_MaxKinFitConstraintVertexPoolSize(size_t locMaxKinFitConstraintVertexPoolSize){dMaxKinFitConstraintVertexPoolSize = locMaxKinFitConstraintVertexPoolSize;}
		inline size_t Get_MaxKinFitConstraintVertexPoolSize(void) const{return dMaxKinFitConstraintVertexPoolSize;}
		inline void Set_MaxKinFitConstraintSpacetimePoolSize(size_t locMaxKinFitConstraintSpacetimePoolSize){dMaxKinFitConstraintSpacetimePoolSize = locMaxKinFitConstraintSpacetimePoolSize;}
		inline size_t Get_MaxKinFitConstraintSpacetimePoolSize(void) const{return dMaxKinFitConstraintSpacetimePoolSize;}
		inline void Set_MaxKinFitConstraintP4PoolSize(size_t locMaxKinFitConstraintP4PoolSize){dMaxKinFitConstraintP4PoolSize = locMaxKinFitConstraintP4PoolSize;}
		inline size_t Get_MaxKinFitConstraintP4PoolSize(void) const{return dMaxKinFitConstraintP4PoolSize;}

		inline void Set_MaxMatrixDSymPoolSize(size_t locMaxMatrixDSymPoolSize){dMaxMatrixDSymPoolSize = locMaxMatrixDSymPoolSize;}
		inline size_t Get_MaxMatrixDSymPoolSize(void) const{return dMaxMatrixDSymPoolSize;}

		//GET FIT INFORMATION
		inline unsigned int Get_NumUnknowns(void) const{return dNumXi;}
		inline unsigned int Get_NumMeasurables(void) const{return dNumEta;}
		inline unsigned int Get_NumConstraintEquations(void) const{return dNumF;}

		//GET FIT QUALITY RESULTS
		inline double Get_ChiSq(void) const{return dChiSq;}
		inline double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}
		inline unsigned int Get_NDF(void) const{return dNDF;}
		inline void Get_Pulls(map<const DKinFitParticle*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;} //key is particle (NULL for rf-bunch), 2nd key is param type

		//GET FIT PARAMETERS
		inline void Get_KinFitParticleIOMap(map<const DKinFitParticle*, const DKinFitParticle*>& locKinFitParticleIOMap) const
		{
			locKinFitParticleIOMap.clear();
			map<const DKinFitParticle*, DKinFitParticle*>::const_iterator locIterator;
			for(locIterator = dKinFitParticleIOMap.begin(); locIterator != dKinFitParticleIOMap.end(); ++locIterator)
				locKinFitParticleIOMap[locIterator->first] = locIterator->second;
		}
		inline const DKinFitParticle* Get_OutputKinFitParticle(const DKinFitParticle* locInputKinFitParticle) const
		{
			map<const DKinFitParticle*, DKinFitParticle*>::const_iterator locIterator = dKinFitParticleIOMap.find(locInputKinFitParticle);
			return ((locIterator != dKinFitParticleIOMap.end()) ? locIterator->second : NULL);
		}
		inline const DKinFitParticle* Get_InputKinFitParticle(const DKinFitParticle* locOutputKinFitParticle) const
		{
			map<const DKinFitParticle*, DKinFitParticle*>::const_iterator locIterator = dKinFitParticleIOMap.begin();
			for(; locIterator != dKinFitParticleIOMap.end(); ++locIterator)
			{
				if(locIterator->second == locOutputKinFitParticle)
					return locIterator->first;
			}
			return NULL;
		}
		inline void Get_KinFitConstraints(deque<const DKinFitConstraint*>& locKinFitConstraints) const
		{
			for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
				locKinFitConstraints.push_back(dKinFitConstraints[loc_i]);
		}
		inline double Get_RFTime(void) const{return dRFTime;}

		//GET UNCERTAINTIES
		inline const TMatrixDSym* Get_VEta(void) {return dVEta;}
		inline const TMatrixDSym* Get_VXi(void) {return dVXi;}
		inline const TMatrixDSym* Get_V(void) {return dV;}

		virtual bool Get_IsBFieldNearBeamline(void) const = 0;

		bool Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixDSym& locCovarianceMatrix) const;
		bool Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixDSym& locCovarianceMatrix, pair<double, double>& locPathLengthPair) const;
		TMatrixDSym* Get_MatrixDSymResource(void);

	protected:
		virtual TVector3 Get_BField(const TVector3& locPosition) const = 0; //must return in units of Tesla!!

	private:

		DKinFitParticle* GetOrCreate_ClonedParticle(const DKinFitParticle* locKinFitParticle);
		void Clone_KinFitParticles(const deque<const DKinFitParticle*>& locKinFitParticles, deque<DKinFitParticle*>& locClonedKinFitParticles);
		DKinFitParticle* Clone_KinFitParticle(const DKinFitParticle* locKinFitParticle);
		TMatrixDSym* Clone_MatrixDSym(const TMatrixDSym* locMatrix);

		DKinFitConstraint_P4* Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint, bool locCloneParticlesFlag);
		DKinFitConstraint_Vertex* Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint, bool locCloneParticlesFlag);
		DKinFitConstraint_Spacetime* Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint, bool locCloneParticlesFlag);

		DKinFitParticle* Get_KinFitParticleResource(void);
		DKinFitConstraint_Vertex* Get_KinFitConstraintVertexResource(void);
		DKinFitConstraint_Spacetime* Get_KinFitConstraintSpacetimeResource(void);
		DKinFitConstraint_P4* Get_KinFitConstraintP4Resource(void);

		bool Calc_dS(void);
		bool Calc_dU(void);

		void Set_MatrixSizes(void);
		void Resize_Matrices(void);
		void Zero_Matrices(void);

		// PREPARE CONSTRAINTS
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
		map<const DKinFitParticle*, DKinFitParticle*> dKinFitParticleIOMap; //key = input, value = cloned (output)

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

		size_t dMaxKinFitParticlePoolSize;
		size_t dMaxKinFitConstraintVertexPoolSize;
		size_t dMaxKinFitConstraintSpacetimePoolSize;
		size_t dMaxKinFitConstraintP4PoolSize;
		size_t dMaxMatrixDSymPoolSize;

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

#endif // _DKinFitter_

