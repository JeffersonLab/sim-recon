#ifndef _DKinFitUtils_
#define _DKinFitUtils_

#include <deque>
#include <set>
#include <map>
#include <algorithm>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"

#include "DKinFitter.h"
#include "DKinFitChain.h"
#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"
#include "DKinFitConstraint_Mass.h"
#include "DKinFitConstraint_P4.h"
#include "DKinFitConstraint_Vertex.h"
#include "DKinFitConstraint_Spacetime.h"

using namespace std;

class DKinFitUtils //purely virtual: cannot directly instantiate class, can only inherit from it
{
	friend class DKinFitter;

	public:

		/***************************************************************** INITIALIZE ***************************************************************/

		//STRUCTORS
		DKinFitUtils(void);
		virtual ~DKinFitUtils(void);

		//RESET: IF YOU OVERRIDE THESE IN THE DERIVED CLASS, BE SURE TO CALL THE BASE CLASS FUNCTIONS!
		virtual void Reset_NewEvent(void);
		virtual void Reset_NewFit(void){};

		void Preallocate_MatrixMemory(void);

		/************************************************************ CONTROL AND MAPPING ***********************************************************/

		//GET CONTROL
		bool Get_LinkVerticesFlag(void) const{return dLinkVerticesFlag;}
		bool Get_DebugLevel(void) const{return dDebugLevel;}

		//SET CONTROL
		void Set_LinkVerticesFlag(bool locLinkVerticesFlag){dLinkVerticesFlag = locLinkVerticesFlag;}
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}

		//GET INPUT FROM OUTPUT
		DKinFitParticle* Get_InputKinFitParticle(DKinFitParticle* locKinFitParticle) const;

		/************************************************************** CREATE PARTICLES ************************************************************/

		//If multiple constraints, it is EXTREMELY CRITICAL that only one DKinFitParticle be created per particle, so that the particles are correctly linked across constraints!!
		DKinFitParticle* Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		DKinFitParticle* Make_TargetParticle(int locPID, int locCharge, double locMass);

		DKinFitParticle* Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const TMatrixDSym* locCovarianceMatrix);
		DKinFitParticle* Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const TMatrixDSym* locCovarianceMatrix);

		DKinFitParticle* Make_MissingParticle(int locPID, int locCharge, double locMass);
		DKinFitParticle* Make_DecayingParticle(int locPID, int locCharge, double locMass, const set<DKinFitParticle*>& locFromInitialState, const set<DKinFitParticle*>& locFromFinalState);

		/********************************************************** SETUP VERTEX CONSTRAINTS ********************************************************/

		//Does not include initial guesses!
		deque<DKinFitConstraint_Vertex*> Create_VertexConstraints(const DKinFitChain* locKinFitChain, bool locSpacetimeFitFlag);

		/************************************************************* CREATE CONSTRAINTS ***********************************************************/

		DKinFitConstraint_Mass* Make_MassConstraint(DKinFitParticle* locDecayingParticle);
		DKinFitConstraint_P4* Make_P4Constraint(const set<DKinFitParticle*>& locInitialParticles, const set<DKinFitParticle*>& locFinalParticles);
		DKinFitConstraint_Vertex* Make_VertexConstraint(const set<DKinFitParticle*>& locFullConstrainParticles, const set<DKinFitParticle*>& locNoConstrainParticles, TVector3 locVertexGuess = TVector3());
		DKinFitConstraint_Spacetime* Make_SpacetimeConstraint(const set<DKinFitParticle*>& locFullConstrainParticles, const set<DKinFitParticle*>& locOnlyConstrainTimeParticles, 
			const set<DKinFitParticle*>& locNoConstrainParticles, TLorentzVector locSpacetimeGuess = TLorentzVector());

		virtual bool Validate_Constraints(const set<DKinFitConstraint*>& locKinFitConstraints) const; //empty, can override

		/*********************************************************** CALCULATION ROUTINES ***********************************************************/

		//if input flag is true: return the value of the p4 at spot defined by locKinFitParticle->Get_Position() //else at the common vertex
			//useful for setting the momentum: locKinFitParticle->Set_Momentum()
		TLorentzVector Calc_DecayingP4_ByPosition(const DKinFitParticle* locKinFitParticle, bool locAtPositionFlag, bool locDontPropagateAtAllFlag = false) const;

		//if input flag is true: return the value of the p4 at the vertex where the p3-deriving particles are at
			//useful for doing mass constraints
		TLorentzVector Calc_DecayingP4_ByP3Derived(const DKinFitParticle* locKinFitParticle, bool locAtP3DerivedFlag, bool locDontPropagateAtAllFlag = false) const;

		//if input flag is true: return the value of the p4 at the production vertex //else return it at the decay vertex
		TLorentzVector Calc_DecayingP4_ByVertex(const DKinFitParticle* locKinFitParticle, bool locAtProductionVertexFlag, bool locDontPropagateAtAllFlag = false) const;

		bool Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, TMatrixDSym& locCovarianceMatrix) const;

		/********************************************************* BUILD OUTPUT DKINFITCHAIN ********************************************************/

		const DKinFitChain* Build_OutputKinFitChain(const DKinFitChain* locInputKinFitChain, set<DKinFitParticle*>& locKinFitOutputParticles);

		/************************************************************ RESOURCE POOL SIZES ***********************************************************/

		//GET CURRENT POOL SIZES
		size_t Get_KinFitParticlePoolSize(void) const{return dKinFitParticlePool_All.size();};
		size_t Get_KinFitConstraintVertexPoolSize(void) const{return dKinFitConstraintVertexPool_All.size();};
		size_t Get_KinFitConstraintSpacetimePoolSize(void) const{return dKinFitConstraintSpacetimePool_All.size();};
		size_t Get_KinFitConstraintP4PoolSize(void) const{return dKinFitConstraintP4Pool_All.size();};
		size_t Get_KinFitConstraintMassPoolSize(void) const{return dKinFitConstraintMassPool_All.size();};
		size_t Get_KinFitChainPoolSize(void) const{return dKinFitChainPool_All.size();};
		size_t Get_KinFitChainStepPoolSize(void) const{return dKinFitChainStepPool_All.size();};
		size_t Get_MatrixDSymPoolSize(void) const{return dMatrixDSymPool_All.size();};
		size_t Get_LargeMatrixDSymPoolSize(void) const{return dLargeMatrixDSymPool_All.size();};

		//GET MAX POOL SIZES
		size_t Get_MaxKinFitParticlePoolSize(void) const{return dMaxKinFitParticlePoolSize;}
		size_t Get_MaxKinFitConstraintVertexPoolSize(void) const{return dMaxKinFitConstraintVertexPoolSize;}
		size_t Get_MaxKinFitConstraintSpacetimePoolSize(void) const{return dMaxKinFitConstraintSpacetimePoolSize;}
		size_t Get_MaxKinFitConstraintP4PoolSize(void) const{return dMaxKinFitConstraintP4PoolSize;}
		size_t Get_MaxKinFitConstraintMassPoolSize(void) const{return dMaxKinFitConstraintMassPoolSize;}
		size_t Get_MaxKinFitChainPoolSize(void) const{return dMaxKinFitChainPoolSize;}
		size_t Get_MaxKinFitChainStepPoolSize(void) const{return dMaxKinFitChainStepPoolSize;}
		size_t Get_MaxMatrixDSymPoolSize(void) const{return dMaxMatrixDSymPoolSize;}
		size_t Get_MaxLargeMatrixDSymPoolSize(void) const{return dMaxLargeMatrixDSymPoolSize;}

		//SET MAX POOL SIZES
		void Set_MaxKinFitParticlePoolSize(size_t locMaxKinFitParticlePoolSize){dMaxKinFitParticlePoolSize = locMaxKinFitParticlePoolSize;}
		void Set_MaxKinFitConstraintVertexPoolSize(size_t locMaxKinFitConstraintVertexPoolSize){dMaxKinFitConstraintVertexPoolSize = locMaxKinFitConstraintVertexPoolSize;}
		void Set_MaxKinFitConstraintSpacetimePoolSize(size_t locMaxKinFitConstraintSpacetimePoolSize){dMaxKinFitConstraintSpacetimePoolSize = locMaxKinFitConstraintSpacetimePoolSize;}
		void Set_MaxKinFitConstraintP4PoolSize(size_t locMaxKinFitConstraintP4PoolSize){dMaxKinFitConstraintP4PoolSize = locMaxKinFitConstraintP4PoolSize;}
		void Set_MaxKinFitConstraintMassPoolSize(size_t locMaxKinFitConstraintMassPoolSize){dMaxKinFitConstraintMassPoolSize = locMaxKinFitConstraintMassPoolSize;}
		void Set_MaxKinFitChainPoolSize(size_t locMaxKinFitChainPoolSize){dMaxKinFitChainPoolSize = locMaxKinFitChainPoolSize;}
		void Set_MaxKinFitChainStepPoolSize(size_t locMaxKinFitChainStepPoolSize){dMaxKinFitChainStepPoolSize = locMaxKinFitChainStepPoolSize;}
		void Set_MaxMatrixDSymPoolSize(size_t locMaxMatrixDSymPoolSize){dMaxMatrixDSymPoolSize = locMaxMatrixDSymPoolSize;}
		void Set_MaxLargeMatrixDSymPoolSize(size_t locMaxLargeMatrixDSymPoolSize){dMaxLargeMatrixDSymPoolSize = locMaxLargeMatrixDSymPoolSize;}

		/********************************************************** END RESOURCE POOL SIZES *********************************************************/

	protected:

		/************************************************************* ABSTRACT FUNCTIONS ***********************************************************/

		//MUST DEFINE IN A DERIVED CLASS
		//FOR SETUP:
		virtual bool Get_IncludeBeamlineInVertexFitFlag(void) const = 0;
		virtual bool Get_IsDetachedVertex(int locPDG_PID) const = 0;

		//FOR FITS:
		virtual TVector3 Get_BField(const TVector3& locPosition) const = 0; //must return in units of Tesla!!
		virtual bool Get_IsBFieldNearBeamline(void) const = 0;

		/*************************************************************** GET RESOURCES **************************************************************/

		DKinFitChain* Get_KinFitChainResource(void);
		DKinFitChainStep* Get_KinFitChainStepResource(void);
		TMatrixDSym* Get_MatrixDSymResource(void);

		/************************************************************** CLONE CONSTRAINTS ***********************************************************/

		//if need to modify a constraint without disrupting the original: note that particles aren't cloned!
		DKinFitConstraint_P4* Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint);
		DKinFitConstraint_Mass* Clone_KinFitConstraint_Mass(const DKinFitConstraint_Mass* locConstraint);
		DKinFitConstraint_Vertex* Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint);
		DKinFitConstraint_Spacetime* Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint);

		/************************************************************** PROTECTED MEMBERS ***********************************************************/

		DKinFitter* dKinFitter; //is set by DKinFitter constructor!
		bool dLinkVerticesFlag;
		int dDebugLevel;

	private:

		/*************************************************************** GET RESOURCES **************************************************************/

		DKinFitParticle* Get_KinFitParticleResource(void);
		DKinFitConstraint_Vertex* Get_KinFitConstraintVertexResource(void);
		DKinFitConstraint_Spacetime* Get_KinFitConstraintSpacetimeResource(void);
		DKinFitConstraint_P4* Get_KinFitConstraintP4Resource(void);
		DKinFitConstraint_Mass* Get_KinFitConstraintMassResource(void);
		TMatrixDSym* Get_LargeMatrixDSymResource(void);

		/************************************************************** CLONE RESOURCES *************************************************************/

		DKinFitParticle* Clone_KinFitParticle(DKinFitParticle* locKinFitParticle);
		TMatrixDSym* Clone_MatrixDSym(const TMatrixDSym* locMatrix);

		set<DKinFitParticle*> Build_CloneParticleSet(const set<DKinFitParticle*>& locInputParticles, const map<DKinFitParticle*, DKinFitParticle*>& locCloneIOMap) const;
		set<DKinFitConstraint*> Clone_ParticlesAndConstraints(const set<DKinFitConstraint*>& locInputConstraints);

		/********************************************************** SETUP VERTEX CONSTRAINTS ********************************************************/

		//Does not include initial guesses!
		deque<set<DKinFitParticle*> > Setup_VertexConstraints(const DKinFitChain* locKinFitChain);
		deque<DKinFitConstraint_Vertex*> Create_VertexConstraints(const deque<set<DKinFitParticle*> >& locAllVertices, bool locSpacetimeFitFlag);
		void Setup_VertexConstraint(const DKinFitChain* locKinFitChain, size_t locStepIndex, set<DKinFitParticle*>& locVertexParticles, set<size_t>& locIncludedStepIndices);
		void Group_VertexParticles(const set<DKinFitParticle*>& locVertexParticles, set<DKinFitParticle*>& locFullConstrainParticles, set<DKinFitParticle*>& locDecayingParticles, set<DKinFitParticle*>& locOnlyConstrainTimeParticles, set<DKinFitParticle*>& locNoConstrainParticles) const;

		/*********************************************************** CALCULATION ROUTINES ***********************************************************/

		//Don't call directly: Rather, call the public wrappers (simpler)
		TLorentzVector Calc_DecayingP4(const DKinFitParticle* locKinFitParticle, bool locIsConstrainedParticle, double locStateSignMultiplier, bool locDontPropagateAtAllFlag = false) const;

		bool Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixDSym& locCovarianceMatrix, pair<double, double>& locPathLengthPair) const;
		void Calc_DecayingParticleJacobian(const DKinFitParticle* locKinFitParticle, bool locDontPropagateDecayingP3Flag, double locStateSignMultiplier, int locNumEta, TMatrixD& locJacobian) const;

		/*************************************************************** CLONE MAPPING **************************************************************/

		//Cannot map input -> output: many outputs for a given input (same particle used in multiple kinfits)
		map<DKinFitParticle*, DKinFitParticle*> dParticleMap_OutputToInput;

		/************************************************************* CONSTRAINT MAPPING ***********************************************************/

		class DSpacetimeParticles //used for map only
		{
			public:
				DSpacetimeParticles(const set<DKinFitParticle*>& locFullConstrainParticles, const set<DKinFitParticle*>& locOnlyConstrainTimeParticles, const set<DKinFitParticle*>& locNoConstrainParticles) : 
				dFullConstrainParticles(locFullConstrainParticles), dOnlyConstrainTimeParticles(locOnlyConstrainTimeParticles), dNoConstrainParticles(locNoConstrainParticles) {}

				bool operator<(const DSpacetimeParticles& locSpacetimeParticles) const;

				set<DKinFitParticle*> dFullConstrainParticles;
				set<DKinFitParticle*> dOnlyConstrainTimeParticles;
				set<DKinFitParticle*> dNoConstrainParticles;
		};

		map<DKinFitParticle*, DKinFitConstraint_Mass*> dMassConstraintMap;
		map<pair<set<DKinFitParticle*>, set<DKinFitParticle*> >, DKinFitConstraint_P4*> dP4ConstraintMap; //pair: initial/final state
		map<pair<set<DKinFitParticle*>, set<DKinFitParticle*> >, DKinFitConstraint_Vertex*> dVertexConstraintMap; //pair: full/no constrain
		map<DSpacetimeParticles, DKinFitConstraint_Spacetime*> dSpacetimeConstraintMap;

		/************************************************************** RESOURCE POOLS **************************************************************/

		//MAX POOL SIZES
		size_t dMaxKinFitParticlePoolSize;
		size_t dMaxKinFitConstraintVertexPoolSize;
		size_t dMaxKinFitConstraintSpacetimePoolSize;
		size_t dMaxKinFitConstraintP4PoolSize;
		size_t dMaxKinFitConstraintMassPoolSize;
		size_t dMaxKinFitChainPoolSize;
		size_t dMaxKinFitChainStepPoolSize;
		size_t dMaxMatrixDSymPoolSize;
		size_t dMaxLargeMatrixDSymPoolSize;

		//PARTICLES
		deque<DKinFitParticle*> dKinFitParticlePool_All;
		deque<DKinFitParticle*> dKinFitParticlePool_Available;

		//CONSTRAINTS
		deque<DKinFitConstraint_Vertex*> dKinFitConstraintVertexPool_All;
		deque<DKinFitConstraint_Vertex*> dKinFitConstraintVertexPool_Available;

		deque<DKinFitConstraint_Spacetime*> dKinFitConstraintSpacetimePool_All;
		deque<DKinFitConstraint_Spacetime*> dKinFitConstraintSpacetimePool_Available;

		deque<DKinFitConstraint_P4*> dKinFitConstraintP4Pool_All;
		deque<DKinFitConstraint_P4*> dKinFitConstraintP4Pool_Available;

		deque<DKinFitConstraint_Mass*> dKinFitConstraintMassPool_All;
		deque<DKinFitConstraint_Mass*> dKinFitConstraintMassPool_Available;

		//DKINFITCHAIN
		deque<DKinFitChain*> dKinFitChainPool_All;
		deque<DKinFitChain*> dKinFitChainPool_Available;

		deque<DKinFitChainStep*> dKinFitChainStepPool_All;
		deque<DKinFitChainStep*> dKinFitChainStepPool_Available;

		//MATRICES
		deque<TMatrixDSym*> dMatrixDSymPool_All;
		deque<TMatrixDSym*> dMatrixDSymPool_Available;

		deque<TMatrixDSym*> dLargeMatrixDSymPool_All;
		deque<TMatrixDSym*> dLargeMatrixDSymPool_Available;
};

inline DKinFitParticle* DKinFitUtils::Get_InputKinFitParticle(DKinFitParticle* locOutputKinFitParticle) const
{
	map<DKinFitParticle*, DKinFitParticle*>::const_iterator locIterator = dParticleMap_OutputToInput.find(locOutputKinFitParticle);
	return (locIterator != dParticleMap_OutputToInput.end() ? locIterator->second : NULL);
}

inline bool DKinFitUtils::DSpacetimeParticles::operator<(const DKinFitUtils::DSpacetimeParticles& locSpacetimeParticles) const
{
	if(dFullConstrainParticles < locSpacetimeParticles.dFullConstrainParticles)
		return true;
	else if(dFullConstrainParticles > locSpacetimeParticles.dFullConstrainParticles)
		return false;

	if(dOnlyConstrainTimeParticles < locSpacetimeParticles.dOnlyConstrainTimeParticles)
		return true;
	else if(dOnlyConstrainTimeParticles > locSpacetimeParticles.dOnlyConstrainTimeParticles)
		return false;

	return (dNoConstrainParticles < locSpacetimeParticles.dNoConstrainParticles);
}

#endif // _DKinFitUtils_
