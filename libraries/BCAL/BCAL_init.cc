// $Id$

#include "DEventLoop.h"
#include "DFactory_DBCALHit.h"

derror_t BCAL_init(DEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DFactory_DBCALHit());

	return NOERROR;
}
