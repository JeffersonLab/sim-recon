// $Id$

#include "DEventLoop.h"
#include "DFactory_DFDCHit.h"

derror_t FDC_init(DEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new DFactory_DFDCHit());

	return NOERROR;
}
