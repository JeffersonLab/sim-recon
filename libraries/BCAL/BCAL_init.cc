// $Id$

#include "DEventLoop.h"
#include "DFactory_DBCALHit.h"
#include "DFactory_DHDDMBCALHit.h"
#include "DFactory_DBCALMCResponse.h"
#include "DFactory_DBCALGeometry.h"
#include "DFactory_DBCALShower.h"

derror_t BCAL_init(DEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DFactory_DBCALHit());
	loop->AddFactory(new DFactory_DHDDMBCALHit());
	loop->AddFactory(new DFactory_DBCALMCResponse());
	loop->AddFactory(new DFactory_DBCALGeometry());
	loop->AddFactory(new DFactory_DBCALShower());

	return NOERROR;
}
