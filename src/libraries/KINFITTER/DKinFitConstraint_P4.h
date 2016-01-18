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

		int Get_FIndex(void) const{return dFIndex;}
		set<DKinFitParticle*> Get_InitialParticles(void) const{return dInitialParticles;};
		set<DKinFitParticle*> Get_FinalParticles(void) const{return dFinalParticles;};

		DKinFitParticle* Get_MissingParticle(void) const; //NULL if none
		DKinFitParticle* Get_OpenEndedDecayingParticle(void) const; //NULL if none
		DKinFitParticle* Get_DefinedParticle() const; //missing or open-ended decaying particle

		set<DKinFitParticle*> Get_AllParticles(void) const;
		void Print_ConstraintInfo(void) const;

	private:

		DKinFitConstraint_P4(void);
		~DKinFitConstraint_P4(void){}

		void Reset(void);
		void Set_FIndex(int locFIndex){dFIndex = locFIndex;}

		void Set_InitialParticles(const set<DKinFitParticle*>& locInitialParticles){dInitialParticles = locInitialParticles;}
		void Set_FinalParticles(const set<DKinFitParticle*>& locFinalParticles){dFinalParticles = locFinalParticles;}

		int dFIndex; //starting row index of the equation(s) corresponding to these constraints in the dF matrix term
		set<DKinFitParticle*> dInitialParticles;
		set<DKinFitParticle*> dFinalParticles;

		TVector3 dInitP3Guess; //initial guess for missing or open-ended-decaying particle. ignored if not present
};

DKinFitConstraint_P4::DKinFitConstraint_P4(void)
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

inline set<DKinFitParticle*> DKinFitConstraint_P4::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	set_union(dInitialParticles.begin(), dInitialParticles.end(), dFinalParticles.begin(), dFinalParticles.end(), inserter(locAllParticles, locAllParticles.begin()));
	return locAllParticles;
}

inline DKinFitParticle* DKinFitConstraint_P4::Get_MissingParticle(void) const
{
	set<DKinFitParticle*> locAllParticles = Get_AllParticles();
	set<DKinFitParticle*>::const_iterator locIterator = locAllParticles.begin();
	for(; locIterator != locAllParticles.end(); ++locIterator)
	{
		DKinFitParticle* locKinFitParticle = *locIterator;
		if(locKinFitParticle->Get_KinFitParticleType() == d_MissingParticle)
			return locKinFitParticle;
	}
	return NULL;
}

inline DKinFitParticle* DKinFitConstraint_P4::Get_OpenEndedDecayingParticle(void) const
{
	//look for decaying particles
	set<DKinFitParticle*> locAllParticles = Get_AllParticles();
	set<DKinFitParticle*>::iterator locIterator = locAllParticles.begin();
	for(; locIterator != locAllParticles.end(); ++locIterator)
	{
		DKinFitParticle* locKinFitParticle = *locIterator;
		if(locKinFitParticle->Get_KinFitParticleType() != d_DecayingParticle)
			continue;

		//see if any of the defined-from particles match those in the constraint
		set<DKinFitParticle*> locFromFinalState = locKinFitParticle->Get_FromFinalState();
		set<DKinFitParticle*> locMatchingParticles;
		set_intersection(locFromFinalState.begin(), locFromFinalState.end(), locAllParticles.begin(), 
			locAllParticles.end(), inserter(locMatchingParticles, locMatchingParticles.begin()));
		if(!locMatchingParticles.empty())
			return locKinFitParticle; //open-ended decaying particle
	}
	return NULL;
}

inline DKinFitParticle* DKinFitConstraint_P4::Get_DefinedParticle() const
{
	DKinFitParticle* locKinFitParticle = Get_MissingParticle();
	if(locKinFitParticle != NULL)
		return locKinFitParticle;
	return Get_OpenEndedDecayingParticle();
}

inline void DKinFitConstraint_P4::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_P4: Initial-state particle PID's, q's, masses: " << endl;
	set<DKinFitParticle*>::const_iterator locIterator = dInitialParticles.begin();
	for(; locIterator != dInitialParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << (*locIterator)->Get_Charge() << ", " << (*locIterator)->Get_Mass() << endl;

	cout << "DKinFitConstraint_P4: Final-state particle PID's, q's, masses: " << endl;
	for(locIterator = dFinalParticles.begin(); locIterator != dFinalParticles.end(); ++locIterator)
		cout << (*locIterator)->Get_PID() << ", " << (*locIterator)->Get_Charge() << ", " << (*locIterator)->Get_Mass() << endl;
}

#endif // _DKinFitConstraint_P4_

