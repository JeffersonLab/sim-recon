#ifndef _DKinFitUtils_
#define _DKinFitUtils_

#include <deque>
#include <set>
#include <map>
#include <algorithm>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixFSym.h"

#include "DResourcePool.h"

#include "DKinFitter.h"
#include "DKinFitChain.h"
#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"
#include "DKinFitConstraint_Mass.h"
#include "DKinFitConstraint_P4.h"
#include "DKinFitConstraint_Vertex.h"
#include "DKinFitConstraint_Spacetime.h"

using namespace std;

class DKinFitUtils //contains pure-virtual functions: cannot directly instantiate class, can only inherit from it
{
	friend class DKinFitter;

	public:

		/***************************************************************** INITIALIZE ***************************************************************/

		//STRUCTORS
		DKinFitUtils(void);
		virtual ~DKinFitUtils(void){};

		//RESET: IF YOU OVERRIDE THESE IN THE DERIVED CLASS, BE SURE TO CALL THE BASE CLASS FUNCTIONS!
		virtual void Reset_NewEvent(void);
		virtual void Reset_NewFit(void){};

		/************************************************************ CONTROL AND MAPPING ***********************************************************/

		//GET CONTROL
		bool Get_LinkVerticesFlag(void) const{return dLinkVerticesFlag;}
		bool Get_DebugLevel(void) const{return dDebugLevel;}
		bool Get_UpdateCovarianceMatricesFlag(void) const{return dUpdateCovarianceMatricesFlag;}

		//SET CONTROL
		void Set_LinkVerticesFlag(bool locLinkVerticesFlag){dLinkVerticesFlag = locLinkVerticesFlag;}
		void Set_DebugLevel(int locDebugLevel){dDebugLevel = locDebugLevel;}
		void Set_UpdateCovarianceMatricesFlag(bool locUpdateCovarianceMatricesFlag){dUpdateCovarianceMatricesFlag = locUpdateCovarianceMatricesFlag;}

		//GET INPUT FROM OUTPUT
		shared_ptr<DKinFitParticle> Get_InputKinFitParticle(const shared_ptr<DKinFitParticle>& locKinFitParticle) const;

		/************************************************************** CREATE PARTICLES ************************************************************/

		//If multiple constraints, it is EXTREMELY CRITICAL that only one DKinFitParticle be created per particle, so that the particles are correctly linked across constraints!!
		shared_ptr<DKinFitParticle> Make_BeamParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix);
		shared_ptr<DKinFitParticle> Make_TargetParticle(int locPID, int locCharge, double locMass);

		//locPathLength is from timing detector to vertex (for updating the time)
		shared_ptr<DKinFitParticle> Make_DetectedParticle(int locPID, int locCharge, double locMass, TLorentzVector locSpacetimeVertex, TVector3 locMomentum, double locPathLength, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix);
		shared_ptr<DKinFitParticle> Make_DetectedShower(int locPID, double locMass, TLorentzVector locSpacetimeVertex, double locShowerEnergy, const shared_ptr<const TMatrixFSym>& locCovarianceMatrix);

		shared_ptr<DKinFitParticle> Make_MissingParticle(int locPID, int locCharge, double locMass);
		shared_ptr<DKinFitParticle> Make_DecayingParticle(int locPID, int locCharge, double locMass, const set<shared_ptr<DKinFitParticle>>& locFromInitialState, const set<shared_ptr<DKinFitParticle>>& locFromFinalState);

		/************************************************************* CREATE CONSTRAINTS ***********************************************************/

		shared_ptr<DKinFitConstraint_Mass> Make_MassConstraint(const shared_ptr<DKinFitParticle>& locDecayingParticle);
		shared_ptr<DKinFitConstraint_P4> Make_P4Constraint(const set<shared_ptr<DKinFitParticle>>& locInitialParticles, const set<shared_ptr<DKinFitParticle>>& locFinalParticles);
		shared_ptr<DKinFitConstraint_Vertex> Make_VertexConstraint(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles, const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles, TVector3 locVertexGuess = TVector3());
		shared_ptr<DKinFitConstraint_Spacetime> Make_SpacetimeConstraint(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles, const set<shared_ptr<DKinFitParticle>>& locOnlyConstrainTimeParticles,
			const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles, TLorentzVector locSpacetimeGuess = TLorentzVector());

		virtual bool Validate_Constraints(const set<shared_ptr<DKinFitConstraint>>& locKinFitConstraints) const; //empty, can override

		/*********************************************************** CALCULATION ROUTINES ***********************************************************/

		//if input flag is true: return the value of the p4 at spot defined by locKinFitParticle->Get_Position() //else at the common vertex
			//useful for setting the momentum: locKinFitParticle->Set_Momentum()
		TLorentzVector Calc_DecayingP4_ByPosition(const DKinFitParticle* locKinFitParticle, bool locAtPositionFlag, bool locDontPropagateAtAllFlag = false) const;

		//if input flag is true: return the value of the p4 at the vertex where the p3-deriving particles are at
			//useful for doing mass constraints
		TLorentzVector Calc_DecayingP4_ByP3Derived(const DKinFitParticle* locKinFitParticle, bool locAtP3DerivedFlag, bool locDontPropagateAtAllFlag = false) const;

		//if input flag is true: return the value of the p4 at the production vertex //else return it at the decay vertex
		TLorentzVector Calc_DecayingP4_ByVertex(const DKinFitParticle* locKinFitParticle, bool locAtProductionVertexFlag, bool locDontPropagateAtAllFlag = false) const;

		bool Propagate_TrackInfoToCommonVertex(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, TVector3& locMomentum, TLorentzVector& locSpacetimeVertex, pair<double, double>& locPathLengthPair, pair<double, double>& locRestFrameLifetimePair, TMatrixFSym* locCovarianceMatrix) const;

		/********************************************************** DKINFITCHAIN RESOURCES **********************************************************/

		//Build output chain
		shared_ptr<const DKinFitChain> Build_OutputKinFitChain(const shared_ptr<const DKinFitChain>& locInputKinFitChain, set<shared_ptr<DKinFitParticle>>& locKinFitOutputParticles);

	protected:

		void Print_Matrix(const TMatrixD& locMatrix) const;
		void Print_Matrix(const TMatrixF& locMatrix) const;

		/************************************************************* ABSTRACT FUNCTIONS ***********************************************************/

		//MUST DEFINE IN A DERIVED CLASS
		//FOR SETUP:
		virtual bool Get_IncludeBeamlineInVertexFitFlag(void) const = 0;

		//FOR FITS:
		virtual TVector3 Get_BField(const TVector3& locPosition) const = 0; //must return in units of Tesla!!
		virtual bool Get_IsBFieldNearBeamline(void) const = 0;

		/********************************************************* GET AND RECYCLE RESOURCES ********************************************************/

		shared_ptr<TMatrixFSym> Get_SymMatrixResource(unsigned int locNumMatrixRows);

		/************************************************************** CLONE RESOURCES *************************************************************/

		//if need to modify a constraint without disrupting the original: note that particles aren't cloned!
		shared_ptr<DKinFitConstraint_P4> Clone_KinFitConstraint_P4(const DKinFitConstraint_P4* locConstraint);
		shared_ptr<DKinFitConstraint_Mass> Clone_KinFitConstraint_Mass(const DKinFitConstraint_Mass* locConstraint);
		shared_ptr<DKinFitConstraint_Vertex> Clone_KinFitConstraint_Vertex(const DKinFitConstraint_Vertex* locConstraint);
		shared_ptr<DKinFitConstraint_Spacetime> Clone_KinFitConstraint_Spacetime(const DKinFitConstraint_Spacetime* locConstraint);

		shared_ptr<TMatrixFSym> Clone_SymMatrix(const TMatrixFSym* locMatrix); //use sparingly in inherited class (if at all)!!

		/************************************************************* RECYCLE RESOURCES ************************************************************/

		// Do this if you are discarding the results from the previous fit (e.g. fit failed, or used to get a vertex guess)
		// Functions are virtual in case the inheriting class wants to manage the memory differently
		void Recycle_LastFitMemory(set<shared_ptr<DKinFitConstraint>>& locKinFitConstraints);

		/************************************************************** PROTECTED MEMBERS ***********************************************************/

		bool Get_IsDecayingParticleDefinedByProducts(const DKinFitParticle* locKinFitParticle) const;

		DKinFitter* dKinFitter; //is set by DKinFitter constructor!
		bool dLinkVerticesFlag;
		int dDebugLevel;
		bool dUpdateCovarianceMatricesFlag;

		shared_ptr<DResourcePool<DKinFitChainStep>> dResourcePool_KinFitChainStep;
		shared_ptr<DResourcePool<DKinFitChain>> dResourcePool_KinFitChain;

	private:

		/************************************************************** CLONE RESOURCES *************************************************************/

		shared_ptr<DKinFitParticle> Clone_KinFitParticle(const shared_ptr<DKinFitParticle>& locKinFitParticle);
		set<shared_ptr<DKinFitParticle>> Build_CloneParticleSet(const set<shared_ptr<DKinFitParticle>>& locInputParticles, const map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitParticle>>& locCloneIOMap) const;
		set<shared_ptr<DKinFitConstraint>> Clone_ParticlesAndConstraints(const set<shared_ptr<DKinFitConstraint>>& locInputConstraints);

		/*********************************************************** CALCULATION ROUTINES ***********************************************************/

		//Don't call directly: Rather, call the public wrappers (simpler)
		TLorentzVector Calc_DecayingP4(const DKinFitParticle* locKinFitParticle, bool locIsConstrainedParticle, double locStateSignMultiplier, bool locDontPropagateAtAllFlag = false) const;

		bool Calc_PathLength(const DKinFitParticle* locKinFitParticle, const TMatrixDSym* locVXi, const TMatrixFSym* locCovarianceMatrix, pair<double, double>& locPathLengthPair, pair<double, double>& locRestFrameLifetimePair) const;
		void Calc_DecayingParticleJacobian(const DKinFitParticle* locKinFitParticle, bool locDontPropagateDecayingP3Flag, double locStateSignMultiplier, int locNumEta, const map<const DKinFitParticle*, int>& locAdditionalPxParamIndices, TMatrixD& locJacobian) const;

		/*************************************************************** CLONE MAPPING **************************************************************/

		//Cannot map input -> output: many outputs for a given input (same particle used in multiple kinfits)
		map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitParticle>> dParticleMap_OutputToInput;

		/************************************************************* CONSTRAINT MAPPING ***********************************************************/

		class DSpacetimeParticles //used for map only
		{
			public:
				DSpacetimeParticles(const set<shared_ptr<DKinFitParticle>>& locFullConstrainParticles, const set<shared_ptr<DKinFitParticle>>& locOnlyConstrainTimeParticles, const set<shared_ptr<DKinFitParticle>>& locNoConstrainParticles) :
				dFullConstrainParticles(locFullConstrainParticles), dOnlyConstrainTimeParticles(locOnlyConstrainTimeParticles), dNoConstrainParticles(locNoConstrainParticles) {}

				bool operator<(const DSpacetimeParticles& locSpacetimeParticles) const;

				set<shared_ptr<DKinFitParticle>> dFullConstrainParticles;
				set<shared_ptr<DKinFitParticle>> dOnlyConstrainTimeParticles;
				set<shared_ptr<DKinFitParticle>> dNoConstrainParticles;
		};

		//Maps of user-created constraints, with the inputs necessary to create them as the keys.
		//These are used to save memory: If a duplicate set of information is entered, instead of creating a new constraint, just return the original one.
		//At the beginning of the fit, these resources are cloned, so that the originals (these) are not modified by the fit.
		//Thus, they can be reused between fits/combos/reactions.
		map<shared_ptr<DKinFitParticle>, shared_ptr<DKinFitConstraint_Mass>> dMassConstraintMap;
		map<pair<set<shared_ptr<DKinFitParticle>>, set<shared_ptr<DKinFitParticle>> >, shared_ptr<DKinFitConstraint_P4>> dP4ConstraintMap; //pair: initial/final state
		map<pair<set<shared_ptr<DKinFitParticle>>, set<shared_ptr<DKinFitParticle>> >, shared_ptr<DKinFitConstraint_Vertex>> dVertexConstraintMap; //pair: full/no constrain
		map<DSpacetimeParticles, shared_ptr<DKinFitConstraint_Spacetime>> dSpacetimeConstraintMap;

		/************************************************************** RESOURCE POOLS **************************************************************/

		shared_ptr<DResourcePool<DKinFitConstraint_Mass>> dResourcePool_MassConstraint;
		shared_ptr<DResourcePool<DKinFitConstraint_P4>> dResourcePool_P4Constraint;
		shared_ptr<DResourcePool<DKinFitConstraint_Vertex>> dResourcePool_VertexConstraint;
		shared_ptr<DResourcePool<DKinFitConstraint_Spacetime>> dResourcePool_SpacetimeConstraint;

		shared_ptr<DResourcePool<TMatrixFSym>> dResourcePool_TMatrixFSym;
		shared_ptr<DResourcePool<DKinFitParticle>> dResourcePool_KinFitParticle;
};

inline shared_ptr<DKinFitParticle> DKinFitUtils::Get_InputKinFitParticle(const shared_ptr<DKinFitParticle>& locOutputKinFitParticle) const
{
	auto locIterator = dParticleMap_OutputToInput.find(locOutputKinFitParticle);
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
