#ifndef _DKinFitConstraint_
#define _DKinFitConstraint_

#include <memory>
#include <set>

#include "DResettable.h"
#include "DKinFitParticle.h"

using namespace std;

class DKinFitConstraint : public DResettable //purely virtual: cannot directly instantiate class, can only inherit from it
{
	public:
		virtual set<shared_ptr<DKinFitParticle>> Get_AllParticles(void) const = 0;
		virtual void Print_ConstraintInfo(void) const = 0;

	protected:
		virtual ~DKinFitConstraint(void) = 0; //forces abstractness
};

inline DKinFitConstraint::~DKinFitConstraint(void){}

#endif // _DKinFitConstraint_

