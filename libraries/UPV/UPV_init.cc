// $Id$

#include "DEventLoop.h"
#include "DFactory_DUPVHit.h"

derror_t UPV_init(DEventLoop *loop)
{
	/// Create and register UPV data factories
	loop->AddFactory(new DFactory_DUPVHit());

	return NOERROR;
}
