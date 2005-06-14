// $Id$

#include "DEventLoop.h"
#include "DFactory_DTOFHit.h"

derror_t TOF_init(DEventLoop *loop)
{
	/// Create and register TOF data factories
	loop->AddFactory(new DFactory_DTOFHit());

	return NOERROR;
}
