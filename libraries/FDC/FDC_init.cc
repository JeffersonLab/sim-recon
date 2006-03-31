// $Id$

#include "DEventLoop.h"
#include "DFactory_DFDCHit.h"
#include "DFactory_DFDCTruth.h"

derror_t FDC_init(DEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new DFactory_DFDCHit());
	loop->AddFactory(new DFactory_DFDCTruth());
	

	return NOERROR;
}
