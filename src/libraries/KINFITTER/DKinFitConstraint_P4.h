#ifndef _DKinFitConstraints_P4_
#define _DKinFitConstraints_P4_

#include <set>
#include <algorithm>

#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"

using namespace std;

class DKinFitter;
class DKinFitUtils;

class DKinFitConstraint_P4 : public DKinFitConstraint
{
	friend class DKinFitter;
	friend class DKinFitUtils;

	public:

		TVector3 Get_InitP3Guess(void) const{return dInitP3Guess;};
		void Set_InitP3Guess(const TVector3& locInitP3Guess){dInitP3Guess = locInitP3Guess;};

		char Get_FIndex(void) const{return dFIndex;}
		set<shared_ptr<DKinFitParticle>> Get_InitialParticles(void) const{return dInitialParticles;};
		set<shared_ptr<DKinFitParticle>> Get_FinalParticles(void) const{return dFinalParticles;};

		shared_ptr<DKinFitParticle> Get_MissingParticle(void) const; //NULL if none
		shared_ptr<DKinFitParticle> Get_OpenEndedDecayingParticle(void) const; //NULL if none
		shared_ptr<DKinFitParticle> Get_DefinedParticle(void) const; //missing or open-ended decaying particle
		bool Get_IsDefinedParticleInFinalState(void) const; //false if initial state or no defined particle

		set<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const;
		void Print_ConstraintInfo(void) const;

	private:

		DKinFitConstraint_P4(void);
		~DKinFitConstraint_P4(void){}

		void Reset(void);
		void Set_FIndex(char locFIndex){dFIndex = locFIndex;}

		void Set_InitialParticles(const set<shared_ptr<DKinFitParticle>>& locInitialParticles){dInitialParticles = locInitialParticles;}
		void Set_FinalParticles(const set<shared_ptr<DKinFitParticle>>& locFinalParticles){dFinalParticles = locFinalParticles;}

		char dFIndex; //starting row index of the equation(s) corresponding to these constraints in the dF matrix term
		set<shared_ptr<DKinFitParticle>> dInitialParticles;
		set<shared_ptr<DKinFitParticle>> dFinalParticles;

		TVector3 dInitP3Guess; //initial guess for missing or open-ended-decaying particle. ignored if not present
};

inline DKinFitConstraint_P4::DKinFitConstraint_P4(void)
{
	Reset();
}

inline void DKinFitConstraint_P4::Reset(void)
{
	dFIndex = 0;
	dInitP3Guess = TVector3(0.0, 0.0, 0.0);
	dInitialParticles.clear();
	dFinalParticles.clear();
}

inline set<shared_ptr<DKinFitParticle>> DKinFitConstraint_P4::Get_AllParticles(void) const
{
	set<shared_ptr<DKinFitParticle>> locAllParticles;
	set_union(dInitialParticles.begin(), dInitialParticles.end(), dFinalParticles.begin(), dFinalParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

inline shared_ptr<DKinFitParticle> DKinFitConstraint_P4::Get_MissingParticle(void) const
{
	auto locAllParticles = Get_AllParticles();
	for(auto& locParticle : locAllParticles)
	{
		if(locParticle->Get_KinFitParticleType() == d_MissingParticle)
			return locParticle;
	}
	return NULL;
}

inline shared_ptr<DKinFitParticle> DKinFitConstraint_P4::Get_OpenEndedDecayingParticle(void) const
{
	//look for decaying particles
	auto locAllParticles = Get_AllParticles();
	for(auto& locParticle : locAllParticles)
	{
		if(locParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue;

		//see if any of the defined-from particles match those in the constraint
		auto locFromFinalState = locParticle->Get_FromFinalState();
		set<shared_ptr<DKinFitParticle>> locMatchingParticles;
		set_intersection(locFromFinalState.begin(), locFromFinalState.end(), locAllParticles.begin(), locAllParticles.end(), inserter(locMatchingParticles, locMatchingParticles.begin()));
		if(!locMatchingParticles.empty())
			return locParticle; //open-ended decaying particle
	}
	return NULL;
}

inline shared_ptr<DKinFitParticle> DKinFitConstraint_P4::Get_DefinedParticle(void) const
{
	auto locKinFitParticle = Get_MissingParticle();
	if(locKinFitParticle != NULL)
		return locKinFitParticle;
	return Get_OpenEndedDecayingParticle();
}

inline bool DKinFitConstraint_P4::Get_IsDefinedParticleInFinalState(void) const
{
	//false if initial state or no defined particle
	auto locDefinedParticle = Get_DefinedParticle();
	if(locDefinedParticle == NULL)
		return false;

	return (dFinalParticles.find(locDefinedParticle) != dFinalParticles.end());
}

inline void DKinFitConstraint_P4::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_P4: Initial-state particle PID's, pointers: " << endl;
	for(auto& locParticle : dInitialParticles)
		cout << locParticle->Get_PID() << ", " << locParticle << endl;

	cout << "DKinFitConstraint_P4: Final-state particle PID's, pointers: " << endl;
	for(auto& locParticle : dFinalParticles)
		cout << locParticle->Get_PID() << ", " << locParticle << endl;
}

#endif // _DKinFitConstraint_P4_

