// $Id$

#include "DEventLoop.h"
#include "DFactory_DHDDMFDCHit.h"
#include "DFactory_DHDDMFDCTruth.h"

derror_t FDC_init(DEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new DFactory_DHDDMFDCHit());
	loop->AddFactory(new DFactory_DHDDMFDCTruth());
	

	return NOERROR;
}
