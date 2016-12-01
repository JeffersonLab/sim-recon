#ifndef _DKinFitter_
#define _DKinFitter_

#include <algorithm>
#include <utility>
#include <math.h>
#include <iostream>
#include <map>
#include <set>
#include <limits>

#include "TVector3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixFSym.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TDecompLU.h"

#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"
#include "DKinFitConstraint_Mass.h"
#include "DKinFitConstraint_P4.h"
#include "DKinFitConstraint_Vertex.h"
#include "DKinFitConstraint_Spacetime.h"

using namespace std;

enum DKinFitStatus
{
	d_KinFitSuccessful = 0,
	d_KinFitFailedSetup,
	d_KinFitFailedInversion,
	d_KinFitTooManyIterations
};

class DKinFitUtils; //forward declaration

class DKinFitter //purely virtual: cannot directly instantiate class, can only inherit from it
{
	public:

		/****************************************************************** PRIMARY *****************************************************************/

		//CONSTRUCTOR
		DKinFitter(DKinFitUtils* locKinFitUtils);

		//RESET
		void Reset_NewEvent(void);
		void Reset_NewFit(void);

		//ADD CONSTRAINTS
		void Add_Constraint(DKinFitConstraint* locKinFitConstraint){dKinFitConstraints.insert(locKinFitConstraint);}
		void Add_Constraints(const set<DKinFitConstraint*>& locKinFitConstraints);

		//FIT!
		bool Fit_Reaction(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!

		/**************************************************************** FIT CONTROL ***************************************************************/

		//GET CONTROL VARIABLES
		int Get_DebugLevel(void) const{return dDebugLevel;}
		unsigned int Get_MaxNumIterations(void) const{return dMaxNumIterations;}
		double Get_ConvergenceChiSqDiff(void) const{return dConvergenceChiSqDiff;}
		double Get_ConvergenceChiSqDiff_LastResort(void) const{return dConvergenceChiSqDiff_LastResort;}

		//SET CONTROL VARIABLES
		void Set_DebugLevel(int locDebugLevel);
		void Set_MaxNumIterations(unsigned int locMaxNumIterations){dMaxNumIterations = locMaxNumIterations;}
		void Set_ConvergenceChiSqDiff(double locConvergenceChiSqDiff){dConvergenceChiSqDiff = locConvergenceChiSqDiff;}
		void Set_ConvergenceChiSqDiff_LastResort(double locConvergenceChiSqDiff){dConvergenceChiSqDiff_LastResort = locConvergenceChiSqDiff;}

		/************************************************************** GET FIT RESULTS *************************************************************/

		//GET STATUS (Fit_Reaction returns true/false for success/fail, but this provides more details on the failures
		DKinFitStatus Get_KinFitStatus(void) const{return dKinFitStatus;}

		//GET FIT INFORMATION
		unsigned int Get_NumUnknowns(void) const{return dNumXi;}
		unsigned int Get_NumMeasurables(void) const{return dNumEta;}
		unsigned int Get_NumConstraintEquations(void) const{return dNumF;}

		//GET FIT QUALITY RESULTS
		double Get_ChiSq(void) const{return dChiSq;}
		double Get_ConfidenceLevel(void) const{return dConfidenceLevel;}
		unsigned int Get_NDF(void) const{return dNDF;}
		void Get_Pulls(map<DKinFitParticle*, map<DKinFitPullType, double> >& locPulls) const{locPulls = dPulls;} //key is particle, 2nd key is param type

		//GET UNCERTAINTIES
		const TMatrixDSym& Get_VEta(void) {return dVEta;}
		const TMatrixDSym& Get_VXi(void) {return dVXi;}
		const TMatrixDSym& Get_V(void) {return dV;}

		//GET OUTPUT PARTICLES & CONSTRAINTS
		set<DKinFitConstraint*> Get_KinFitConstraints(void) const{return dKinFitConstraints;}
		set<DKinFitParticle*> Get_KinFitParticles(void) const{return dKinFitParticles;}

		/************************************************************ END GET FIT RESULTS ***********************************************************/

	private:

		//PRIVATE DEFAULT CONSTRUCTOR
		DKinFitter(void){}; //Cannot use default constructor. Must construct with DKinFitUtils as argument

		/************************************************************ UTILITY FUNCTIONS *************************************************************/

		void Print_Matrix(const TMatrixD& locMatrix) const;

		template <typename DType> set<DType*> Get_Constraints(void) const;
		template <typename DType> set<DType*> Get_Constraints(const set<DKinFitConstraint*>& locConstraints) const;

		//When we ask this question, we want to know, is the particle's information used in a constraint of the given type
			//For example, if a decaying particle is used for a vertex constraint, the particle's that DEFINE the decaying particle are also included
		template <typename DType> set<DType*> Get_Constraints(DKinFitParticle* locKinFitParticle, bool locOnlyDirectFlag = false) const;
		template <typename DType> bool Get_IsInConstraint(DKinFitParticle* locKinFitParticle, bool locOnlyDirectFlag = false) const;
		template <typename DType> bool Get_IsIndirectlyInConstraint(DKinFitParticle* locKinFitParticle) const;

		bool Get_IsConstrainingVertex(DKinFitParticle* locKinFitParticle) const;
		bool Get_IsTimeConstrained(DKinFitParticle* locKinFitParticle) const;

		/*********************************************************** FIT INITIALIZATION *************************************************************/

		void Prepare_ConstraintsAndParticles(void);
		void Set_MatrixSizes(void);
		void Resize_Matrices(void);
		void Zero_Matrices(void);
		void Fill_InputMatrices(void);

		/************************************************************ CALCULATE MATRICES ************************************************************/

		bool Iterate(void);

		bool Calc_dS(void);
		bool Calc_dU(void);
		void Calc_dVdEta(void);

		void Calc_dF(void);

		void Calc_dF_P4(int locFIndex, const DKinFitParticle* locKinFitParticle, double locStateSignMultiplier);
		void Calc_dF_MassDerivs(size_t locFIndex, const DKinFitParticle* locKinFitParticle, TLorentzVector locXP4, double locStateSignMultiplier, bool locIsConstrainedParticle);

		void Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier);
		void Calc_dF_Vertex_NotDecaying(size_t locFIndex, const DKinFitParticle* locKinFitParticle);
		void Calc_dF_Vertex_Decaying_Accel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier);
		void Calc_dF_Vertex_Decaying_NonAccel(size_t locFIndex, const DKinFitParticle* locKinFitParticle, const DKinFitParticle* locKinFitParticle_DecayingSource, double locStateSignMultiplier);
		void Calc_Vertex_Params(const DKinFitParticle* locKinFitParticle, double& locJ, TVector3& locQ, TVector3& locM, TVector3& locD);
		TVector3 Calc_VertexParams_P4DerivedAtCommonVertex(const DKinFitParticle* locKinFitParticle);

		/************************************************************* UPDATE & FINAL ***************************************************************/

		void Update_ParticleParams(void);
		void Calc_Pulls(void);
		void Set_FinalTrackInfo(void);

		/***************************************************** FIT CONTROL AND UTILITY VARIABLES ****************************************************/

		DKinFitUtils* dKinFitUtils;

		DKinFitStatus dKinFitStatus;
		unsigned int dDebugLevel;

		unsigned int dMaxNumIterations;

		double dConvergenceChiSqDiff;
		double dConvergenceChiSqDiff_LastResort; //if max # iterations hit, use this for final check (sometimes chisq walks (very slightly) forever without any meaningful change in the variables)

		/******************************************************** CONSTRAINTS AND PARTICLES *********************************************************/

		set<DKinFitConstraint*> dKinFitConstraints;
		set<DKinFitParticle*> dKinFitParticles;
		map<DKinFitParticle*, set<DKinFitConstraint*> > dParticleConstraintMap; //these particles are either directly or indirectly in these constraints
		map<DKinFitParticle*, set<DKinFitConstraint*> > dParticleConstraintMap_Direct; //these particles are directly in these constraints
		//indirect: if only present to define the momentum of a decaying particle

		/*********************************************************** FIT MATRIX VARIABLES ***********************************************************/

		unsigned int dNumXi; //num unknowns
		unsigned int dNumEta; //num measurables
		unsigned int dNumF; //num constraint eqs

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

		TMatrixDSym dVXi; //covariance matrix of dXi
		TMatrixDSym dVEta; //covariance matrix of dEta
		TMatrixDSym dV; //full covariance matrix: dVEta at top-left and dVXi at bottom-right (+ the eta, xi covariance)

		/*************************************************************** FIT RESULTS ****************************************************************/

		double dChiSq;
		unsigned int dNDF;
		double dConfidenceLevel;

		//Pulls: 2nd dimension can be E, x, y, z, t for neutral showers; px, py, pz, x, y, z, t for charged or neutral tracks; or a subset of these
		map<DKinFitParticle*, map<DKinFitPullType, double> > dPulls; //key is particle, 2nd key is param type
};

inline void DKinFitter::Add_Constraints(const set<DKinFitConstraint*>& locKinFitConstraints)
{
	dKinFitConstraints.insert(locKinFitConstraints.begin(), locKinFitConstraints.end());
}

template <typename DType> inline set<DType*> DKinFitter::Get_Constraints(void) const
{
	//Get all constraints of a given type
	return Get_Constraints<DType>(dKinFitConstraints);
}

template <typename DType> inline set<DType*> DKinFitter::Get_Constraints(const set<DKinFitConstraint*>& locConstraints) const
{
	//Get all constraints of a given type
	set<DType*> locTypeConstraints;
	set<DKinFitConstraint*>::const_iterator locConstraintIterator = locConstraints.begin();
	for(; locConstraintIterator != locConstraints.end(); ++locConstraintIterator)
	{
		DType* locConstraint = dynamic_cast<DType*>(*locConstraintIterator);
		if(locConstraint != NULL)
			locTypeConstraints.insert(locConstraint);
	}
	return locTypeConstraints;
}

template <typename DType> inline bool DKinFitter::Get_IsInConstraint(DKinFitParticle* locKinFitParticle, bool locOnlyDirectFlag) const
{
	//Return whether a particle is in a constraint of the given type
	return !Get_Constraints<DType>(locKinFitParticle, locOnlyDirectFlag).empty();
}

template <typename DType> inline bool DKinFitter::Get_IsIndirectlyInConstraint(DKinFitParticle* locKinFitParticle) const
{
	//Return whether a particle is indirectly in a constraint of the given type //in it, but not directly
	return (Get_IsInConstraint<DType>(locKinFitParticle, false) && !Get_IsInConstraint<DType>(locKinFitParticle, true));
}

template <typename DType> inline set<DType*> DKinFitter::Get_Constraints(DKinFitParticle* locKinFitParticle, bool locOnlyDirectFlag) const
{
	//Get all constraints of a given type that a given particle is in: either directly or INDIRECTLY
		//If the particle is used to define the p3 of a decaying particle that is used in a given constraint, it is included
	const map<DKinFitParticle*, set<DKinFitConstraint*> >& locConstraintMap = locOnlyDirectFlag ? dParticleConstraintMap_Direct : dParticleConstraintMap;
	map<DKinFitParticle*, set<DKinFitConstraint*> >::const_iterator locMapIterator = locConstraintMap.find(locKinFitParticle);
	if(locMapIterator == locConstraintMap.end())
		return set<DType*>();

	const set<DKinFitConstraint*>& locParticleConstraints = locMapIterator->second;
	return Get_Constraints<DType>(locParticleConstraints);
}

#endif // _DKinFitter_
