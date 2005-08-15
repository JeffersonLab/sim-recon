// $Id$

#include "DEventLoop.h"
#include "DFactory_DTOFHit.h"
#include "DFactory_DTOFGeometry.h"
#include "DFactory_DTOFMCResponse.h"

derror_t TOF_init(DEventLoop *loop)
{
	/// Create and register TOF data factories
	loop->AddFactory(new DFactory_DTOFMCResponse());
	loop->AddFactory(new DFactory_DTOFHit());
	loop->AddFactory(new DFactory_DTOFGeometry());

	return NOERROR;
}
