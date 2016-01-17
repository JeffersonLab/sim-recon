#ifndef _DKinFitConstraints_Mass_
#define _DKinFitConstraints_Mass_

#include <set>

#include "DKinFitParticle.h"
#include "DKinFitConstraint.h"

using namespace std;

class DKinFitter;
class DKinFitUtils;

class DKinFitConstraint_Mass : public DKinFitConstraint
{
	friend class DKinFitter;
	friend class DKinFitUtils;

	public:

		DKinFitParticle* Get_DecayingParticle(void) const{return dDecayingParticle;}
		set<DKinFitParticle*> Get_AllParticles(void) const;

		int Get_FIndex(void) const{return dFIndex;}
		void Print_ConstraintInfo(void) const;

	private:

		DKinFitConstraint_Mass(void);
		~DKinFitConstraint_Mass(void){}

		void Reset(void);
		void Set_FIndex(int locFIndex){dFIndex = locFIndex;}
		void Set_DecayingParticle(DKinFitParticle* locDecayingParticle){dDecayingParticle = locDecayingParticle;}

		DKinFitParticle* dDecayingParticle;
		int dFIndex; //starting row index of the equation(s) corresponding to these constraints in the dF matrix term
};

DKinFitConstraint_Mass::DKinFitConstraint_Mass(void)
{
	Reset();
}

void DKinFitConstraint_Mass::Reset(void)
{
	dFIndex = 0;
	dDecayingParticle = NULL;
}

inline set<DKinFitParticle*> DKinFitConstraint_Mass::Get_AllParticles(void) const
{
	set<DKinFitParticle*> locAllParticles;
	locAllParticles.insert(dDecayingParticle);
	return locAllParticles;
}

inline void DKinFitConstraint_Mass::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_Mass: Decaying particle PID, q, mass: " << endl;
	cout << dDecayingParticle->Get_PID() << ", " << dDecayingParticle->Get_Charge() << ", " << dDecayingParticle->Get_Mass() << endl;
}

#endif // _DKinFitConstraint_Mass_

