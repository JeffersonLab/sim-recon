#ifndef _DKinFitConstraint_
#define _DKinFitConstraint_

#include <set>

#include "DKinFitParticle.h"

using namespace std;

class DKinFitConstraint //purely virtual: cannot directly instantiate class, can only inherit from it
{
	public:
		virtual set<DKinFitParticle*> Get_AllParticles(void) const = 0;
		virtual void Print_ConstraintInfo(void) const = 0;

	protected:
		virtual ~DKinFitConstraint(void) = 0; //forces abstractness
};

#endif // _DKinFitConstraint_

