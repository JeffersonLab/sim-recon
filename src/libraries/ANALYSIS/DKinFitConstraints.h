#ifndef _DKinFitConstraints_
#define _DKinFitConstraints_

#include <deque>
#include <map>
#include <set>
#include <utility>

#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

class DKinFitter;
class DKinFitParticle;

class DKinFitConstraint //purely virtual: cannot directly instantiate class, can only inherit from it
{
	friend class DKinFitter;
	public:
		string Get_ConstraintString(void) const{return dConstraintString;}
		void Set_ConstraintString(string locConstraintString){dConstraintString = locConstraintString;}
		virtual void Get_AllKinFitParticles(set<const DKinFitParticle*>& locKinFitParticles) const = 0;

	protected:
		virtual ~DKinFitConstraint(void) = 0; //forces abstractness

		string dConstraintString; //e.g. "m_{#it{#Lambda}}
};

class DKinFitConstraint_VertexBase : public DKinFitConstraint //purely virtual: cannot directly instantiate class, can only inherit from it
{
	friend class DKinFitter;

	public:
		virtual TVector3 Get_CommonVertex(void) const = 0; //inheriting classes MUST define this method
		int Get_VxParamIndex(void) const{return dVxParamIndex;}
		int Get_FIndex(const DKinFitParticle* locKinFitParticle) const;

		void Replace_DecayingParticle(const DKinFitParticle* locOriginalParticle, const DKinFitParticle* locTreatAsDetectedParticle);

		virtual void Set_VertexGuess(TVector3& locVertex) = 0; //inheriting classes MUST define this method

		deque<const DKinFitParticle*> Get_FullConstrainParticles(void) const{return deque<const DKinFitParticle*>(dFullConstrainParticles.begin(), dFullConstrainParticles.end());}
		deque<const DKinFitParticle*> Get_NoConstrainParticles(void) const{return deque<const DKinFitParticle*>(dNoConstrainParticles.begin(), dNoConstrainParticles.end());}
		deque<pair<const DKinFitParticle*, bool> > Get_DecayingParticles(void) const{return deque<pair<const DKinFitParticle*, bool> >(dDecayingParticles.begin(), dDecayingParticles.end());}

		virtual void Get_AllKinFitParticles(set<const DKinFitParticle*>& locKinFitParticles) const
		{
			locKinFitParticles.clear();
			for(size_t loc_i = 0; loc_i < dFullConstrainParticles.size(); ++loc_i)
				locKinFitParticles.insert(dFullConstrainParticles[loc_i]);
			for(size_t loc_i = 0; loc_i < dNoConstrainParticles.size(); ++loc_i)
				locKinFitParticles.insert(dNoConstrainParticles[loc_i]);
			for(size_t loc_i = 0; loc_i < dDecayingParticles.size(); ++loc_i)
				locKinFitParticles.insert(dDecayingParticles[loc_i].first);

			set<DKinFitParticle*>::iterator locIterator = dDecayingParticlesToAssign.begin();
			for(; locIterator != dDecayingParticlesToAssign.end(); ++locIterator)
				locKinFitParticles.insert((*locIterator));
		}

	protected:

		DKinFitConstraint_VertexBase(void);
		virtual ~DKinFitConstraint_VertexBase(void){}

		virtual void Reset(void);

		void Set_FIndex(const DKinFitParticle* locKinFitParticle, unsigned int locFIndex){dConstraintEquationParticleMap[locKinFitParticle] = locFIndex;}
		void Set_VxParamIndex(int locVxParamIndex){dVxParamIndex = locVxParamIndex;}

		void Set_FullConstrainParticles(deque<DKinFitParticle*>& locFullConstrainParticles){dFullConstrainParticles = locFullConstrainParticles;}
		void Set_NoConstrainParticles(deque<DKinFitParticle*>& locNoConstrainParticles){dNoConstrainParticles = locNoConstrainParticles;}
		void Set_DecayingParticles(deque<pair<DKinFitParticle*, bool> >& locDecayingParticles)
		{
			dDecayingParticles = locDecayingParticles;
			for(size_t loc_i = 0; loc_i < dDecayingParticles.size(); ++loc_i)
				dDecayingParticlesToAssign.insert(dDecayingParticles[loc_i].first);
		}

		void Add_FullConstrainParticle(DKinFitParticle* locKinFitParticle){dFullConstrainParticles.push_back(locKinFitParticle);}
		void Add_NoConstrainParticle(DKinFitParticle* locKinFitParticle){dNoConstrainParticles.push_back(locKinFitParticle);}

		int dVxParamIndex; //location of the Vx uncertainty in the dVXi covariance matrix term (Vy & Vz are the subsequent terms)

		deque<pair<DKinFitParticle*, bool> > dDecayingParticles; //bool is true if vertex is production vertex / particle in final state, false if decay vertex / initial state
		deque<DKinFitParticle*> dFullConstrainParticles; //charged particles, decaying particles, beam particles (not neutral showers!)
		deque<DKinFitParticle*> dNoConstrainParticles; //missing particles & neutral showers //not used to constrain vertex or time, but fit vertex is set for this particle

		set<DKinFitParticle*> dDecayingParticlesToAssign; //decaying particles not yet assigned to either constrain or be constrained ("no constrain")

		map<const DKinFitParticle*, unsigned int> dConstraintEquationParticleMap; //key is particle (NULL for RF time), value is the constraint equation index (must check whether decaying or not to see what f-matrix it's in!)
};

class DKinFitConstraint_Vertex : public DKinFitConstraint_VertexBase
{
	friend class DKinFitter;

	public:
		TVector3 Get_CommonVertex(void) const{return dCommonVertex;}
		void Set_VertexGuess(TVector3& locVertex){dCommonVertex = locVertex;}

		bool Get_DecayingParticleInInitialStateFlag(const DKinFitParticle* locKinFitParticle) const;

	private:
		DKinFitConstraint_Vertex(void);
		~DKinFitConstraint_Vertex(void){}

		void Reset(void);

		void Set_CommonVertex(TVector3 locCommonVertex){dCommonVertex = locCommonVertex;}

		TVector3 dCommonVertex; //for propagating the track to this point
};

class DKinFitConstraint_Spacetime : public DKinFitConstraint_VertexBase
{
	friend class DKinFitter;

	public:
		void Set_TimeGuess(double locTime){dCommonSpacetimeVertex.SetT(locTime);}
		void Set_VertexGuess(TVector3& locVertex){dCommonSpacetimeVertex.SetVect(locVertex);}
		void Set_SpacetimeGuess(TLorentzVector& locSpacetimeVertex){dCommonSpacetimeVertex = locSpacetimeVertex;}

		TLorentzVector Get_CommonSpacetimeVertex(void) const{return dCommonSpacetimeVertex;}
		TVector3 Get_CommonVertex(void) const{return dCommonSpacetimeVertex.Vect();}
		double Get_CommonTime(void) const{return dCommonSpacetimeVertex.T();}

		bool Get_UseRFTimeFlag(void) const{return dUseRFTimeFlag;}
		const DKinFitParticle* Get_BeamParticle(void) const{return dBeamParticle;}
		int Get_TParamIndex(void) const{return dTParamIndex;}

		deque<const DKinFitParticle*> Get_OnlyConstrainTimeParticles(void) const{return deque<const DKinFitParticle*>(dOnlyConstrainTimeParticles.begin(), dOnlyConstrainTimeParticles.end());} //neutral showers

		bool Get_DecayingParticleInInitialStateFlag(const DKinFitParticle* locKinFitParticle) const;

		void Get_AllKinFitParticles(set<const DKinFitParticle*>& locKinFitParticles) const
		{
			DKinFitConstraint_VertexBase::Get_AllKinFitParticles(locKinFitParticles);
			for(size_t loc_i = 0; loc_i < dOnlyConstrainTimeParticles.size(); ++loc_i)
				locKinFitParticles.insert(dOnlyConstrainTimeParticles[loc_i]);
		}

	private:
		DKinFitConstraint_Spacetime(void);
		~DKinFitConstraint_Spacetime(void){}

		void Reset(void);

		void Set_CommonVertex(TVector3 locCommonVertex){dCommonSpacetimeVertex.SetVect(locCommonVertex);}
		void Set_CommonTime(double locCommonTime){dCommonSpacetimeVertex.SetT(locCommonTime);}
		void Set_CommonSpacetimeVertex(TLorentzVector locCommonSpacetimeVertex){dCommonSpacetimeVertex = locCommonSpacetimeVertex;}

		void Set_UseRFTimeFlag(bool locUseRFTimeFlag){dUseRFTimeFlag = locUseRFTimeFlag;}
		void Set_TParamIndex(int locTParamIndex){dTParamIndex = locTParamIndex;}

		void Set_OnlyConstrainTimeParticles(deque<DKinFitParticle*>& locOnlyConstrainTimeParticles){dOnlyConstrainTimeParticles = locOnlyConstrainTimeParticles;}

		deque<DKinFitParticle*> dOnlyConstrainTimeParticles; //neutral showers //not used to constrain vertex, but fit vertex is used for time constraint

		bool dUseRFTimeFlag; //for time constraint
		TLorentzVector dCommonSpacetimeVertex; //for propagating the track to this point
		int dTParamIndex;
		DKinFitParticle* dBeamParticle; //NULL if not included in this constraint
};

class DKinFitConstraint_P4 : public DKinFitConstraint //one will be implemented as a p4 constraint, the rest as invariant mass constraints
{
	friend class DKinFitter;

	public:
		int Get_FIndex(void) const{return dFIndex;}
		void Set_ConstrainInitialParticleMassFlag(bool locConstrainMassFlag){dConstrainMassFlag = locConstrainMassFlag;}
		bool Get_ConstrainInitialParticleMassFlag(void) const{return dConstrainMassFlag;}
		bool Get_ConstrainedParticleIsInInitialStateFlag(void) const{return dConstrainedParticleIsInInitialStateFlag;}
		bool Get_ConstrainMassByInvariantMassFlag(void) const{return dConstrainMassByInvariantMassFlag;}
		bool Get_IsActualP4ConstraintFlag(void) const{return dIsActualP4ConstraintFlag;}

		const DKinFitParticle* Get_ConstrainedP4Particle(void) const{return dConstrainedP4Particle;}

		deque<const DKinFitParticle*> Get_InitialParticles(void) const{return deque<const DKinFitParticle*>(dInitialParticles.begin(), dInitialParticles.end());}
		deque<const DKinFitParticle*> Get_FinalParticles(void) const{return deque<const DKinFitParticle*>(dFinalParticles.begin(), dFinalParticles.end());}

		void Replace_Particle(const DKinFitParticle* locOriginalParticle, bool locInitialStateFlag, const DKinFitParticle* locTreatAsDetectedParticle);

		void Get_AllKinFitParticles(set<const DKinFitParticle*>& locKinFitParticles) const
		{
			locKinFitParticles.clear();
			for(size_t loc_i = 0; loc_i < dInitialParticles.size(); ++loc_i)
				locKinFitParticles.insert(dInitialParticles[loc_i]);
			for(size_t loc_i = 0; loc_i < dFinalParticles.size(); ++loc_i)
				locKinFitParticles.insert(dFinalParticles[loc_i]);
		}

	private:
		DKinFitConstraint_P4(void);
		~DKinFitConstraint_P4(void){}

		void Reset(void);

		void Set_FIndex(int locFIndex){dFIndex = locFIndex;}
		void Set_ConstrainedP4Particle(DKinFitParticle* locConstrainedP4Particle){dConstrainedP4Particle = locConstrainedP4Particle;}
		void Set_InitialParticles(const deque<DKinFitParticle*>& locInitialParticles){dInitialParticles = locInitialParticles;}
		void Set_FinalParticles(const deque<DKinFitParticle*>& locFinalParticles){dFinalParticles = locFinalParticles;}
		void Set_ConstrainedParticleIsInInitialStateFlag(bool locFlag){dConstrainedParticleIsInInitialStateFlag = locFlag;}

		int dFIndex; //starting row index of the equation(s) corresponding to these constraints in the dF matrix term // -1 if a decaying particle mass constraint not applied

		bool dConstrainMassFlag;
		bool dIsActualP4ConstraintFlag; //e.g. initial particle is beam or open-ended decaying particle AND this is not an inclusive-p4 fit
		bool dConstrainedParticleIsInInitialStateFlag;
		bool dConstrainMassByInvariantMassFlag; //true unless missing particle with unknown mass is a decay product (ignored if dConstrainMassFlag is false)
		DKinFitParticle* dConstrainedP4Particle;
		deque<DKinFitParticle*> dInitialParticles;
		deque<DKinFitParticle*> dFinalParticles;
};

#endif // _DKinFitConstraints_

