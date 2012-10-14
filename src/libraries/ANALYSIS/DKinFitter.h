#ifndef _DKinFitter_
#define _DKinFitter_

#include <deque>
#include <utility>
#include <math.h>
#include <iostream>
#include <map>
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
		const DKinFitParticle* Make_DecayingParticle(int locCharge, double locMass);
		const DKinFitParticle* Make_MissingParticle(int locCharge, double locMass);
		const DKinFitParticle* Make_BeamParticle(int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		const DKinFitParticle* Make_TargetParticle(int locCharge, double locMass);
		const DKinFitParticle* Make_DetectedParticle(int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		const DKinFitParticle* Make_DetectedShower(double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix);

		//INPUT RF TIME AND FIT CONSTRAINTS
		void Set_RFTime(double locRFTime, double locRFUncertainty, const DKinFitParticle* locRFMatchedBeamParticle);
		const DKinFitConstraint_Vertex* Add_VertexConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, TVector3& locVertexGuess);
		const DKinFitConstraint_Spacetime* Add_SpacetimeConstraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles, bool locUseRFTimeFlag, TVector3& locVertexGuess, double locCommonTimeGuess);
		const DKinFitConstraint_P4* Add_P4Constraint(const deque<const DKinFitParticle*>& locInitialParticles, const deque<const DKinFitParticle*>& locFinalParticles);

		//FIT!
		virtual bool Fit_Reaction(void); //IF YOU OVERRIDE THIS METHOD IN THE DERIVED CLASS, MAKE SURE YOU CALL THIS BASE CLASS METHOD!!

		//GET/SET CONTROL VARIABLES

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
		inline void Get_KinFitConstraints(deque<const DKinFitConstraint*>& locKinFitConstraints) const
		{
			for(size_t loc_i = 0; loc_i < dKinFitConstraints.size(); ++loc_i)
				locKinFitConstraints.push_back(dKinFitConstraints[loc_i]);
		}
		inline double Get_RFTime(void) const{return dRFTime;}

		//GET UNCERTAINTIES
		inline const TMatrixDSym* Get_VEta(void) {return dVEta;}
		inline const TMatrixDSym* Get_VXi(void) {return dVXi;}

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

		void Set_Particle(DKinFitParticle* locKinFitParticle);

		DKinFitParticle* Get_KinFitParticleResource(void);
		DKinFitConstraint_Vertex* Get_KinFitConstraintVertexResource(void);
		DKinFitConstraint_Spacetime* Get_KinFitConstraintSpacetimeResource(void);
		DKinFitConstraint_P4* Get_KinFitConstraintP4Resource(void);

		bool Calc_dS(void);
		bool Calc_dVXi(void);

		void Set_MatrixSizes(void);
		void Setup_Matrices(void);
		void Zero_Matrices(void);

		bool Resolve_Constraints(void);
		bool Resolve_DecayingParticleSpacetimeLinks(void);
		bool Find_ConstrainableParticles(deque<pair<DKinFitParticle*, DKinFitConstraint_P4*> >& locConstrainableParticles, const deque<const DKinFitParticle*>& locConstrainedParticles);
		void Constrain_Particle(pair<DKinFitParticle*, DKinFitConstraint_P4*>& locConstrainableParticle);

		void Fill_InputMatrices(void);
		void Update_ParticleParams(void);

		void Print_Matrix(const TMatrixD& locMatrix) const;
		void Print_ParticleParams(const DKinFitParticle* locKinFitParticle) const;

		void Calc_dF(void);
		void Calc_dF_P4(size_t locFIndex, const DKinFitParticle* locKinFitParticle, bool locInitialStateFlag, const DKinFitParticle* locConstrainedParticle);
		void Calc_dF_Vertex(size_t locFIndex, const DKinFitParticle* locKinFitParticle);
		void Calc_dF_Time(size_t& locFIndex, const DKinFitParticle* locKinFitParticle, bool locUseRFTimeFlag);
		void Calc_Pulls(void);
		void Set_FinalTrackInfo(void);

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
		TMatrixD dR;

		TMatrixD dF; //constraint equations
		TMatrixD dEpsilon; //dY - dEta
		TMatrixD dLambda; //lagrange multipliers
		TMatrixD dLambda_T;

		TMatrixD dFdEta; //partial derivative of constraint equations wrst the observables
		TMatrixD dFdEta_T;
		TMatrixD dFdXi; //partial derivative of constraint equations wrst the unmeasurable unknowns
		TMatrixD dFdXi_T;

		TMatrixDSym* dVXi; //covariance matrix of dXi
		TMatrixDSym dVXi_Inverse;
		TMatrixDSym* dVEta; //covariance matrix of dEta

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
};

#endif // _DKinFitter_

