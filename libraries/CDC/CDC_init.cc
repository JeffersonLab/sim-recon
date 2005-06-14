// $Id$

#include "DEventLoop.h"
#include "DFactory_DCDCHit.h"

derror_t CDC_init(DEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new DFactory_DCDCHit());

	return NOERROR;
}
