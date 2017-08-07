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

		shared_ptr<DKinFitParticle> Get_DecayingParticle(void) const{return dDecayingParticle;}
		set<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const{return {dDecayingParticle};}

		char Get_FIndex(void) const{return dFIndex;}
		void Print_ConstraintInfo(void) const;

	private:

		DKinFitConstraint_Mass(void);
		~DKinFitConstraint_Mass(void){}

		void Reset(void);
		void Set_FIndex(char locFIndex){dFIndex = locFIndex;}
		void Set_DecayingParticle(const shared_ptr<DKinFitParticle>& locDecayingParticle){dDecayingParticle = locDecayingParticle;}

		shared_ptr<DKinFitParticle> dDecayingParticle;
		char dFIndex; //starting row index of the equation(s) corresponding to these constraints in the dF matrix term
};

inline DKinFitConstraint_Mass::DKinFitConstraint_Mass(void)
{
	Reset();
}

inline void DKinFitConstraint_Mass::Reset(void)
{
	dFIndex = 0;
	dDecayingParticle = NULL;
}

inline void DKinFitConstraint_Mass::Print_ConstraintInfo(void) const
{
	cout << "DKinFitConstraint_Mass: Decaying particle PID, pointer: " << endl;
	cout << dDecayingParticle->Get_PID() << ", " << dDecayingParticle << endl;
}

#endif // _DKinFitConstraint_Mass_

