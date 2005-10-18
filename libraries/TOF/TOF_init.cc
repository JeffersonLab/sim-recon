// $Id$

#include "DEventLoop.h"
#include "DFactory_DTOFHit.h"
#include "DFactory_DTOFGeometry.h"
#include "DFactory_DTOFMCResponse.h"
#include "DFactory_DHDDMTOFTruth.h"
#include "DFactory_DHDDMTOFHit.h"
#include "DFactory_DTOFPoint.h"

derror_t TOF_init(DEventLoop *loop)
{
	/// Create and register TOF data factories
	loop->AddFactory(new DFactory_DTOFMCResponse());
	loop->AddFactory(new DFactory_DTOFHit());
	loop->AddFactory(new DFactory_DTOFGeometry());
	loop->AddFactory(new DFactory_DHDDMTOFTruth());
	loop->AddFactory(new DFactory_DHDDMTOFHit());
	loop->AddFactory(new DFactory_DTOFPoint());

	return NOERROR;
}
