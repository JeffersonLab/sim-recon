// $Id$

#include "DEventLoop.h"
#include "DFactory_DCHERENKOVHit.h"

derror_t CHERENKOV_init(DEventLoop *loop)
{
	/// Create and register CHERENKOV data factories
	loop->AddFactory(new DFactory_DCHERENKOVHit());

	return NOERROR;
}
