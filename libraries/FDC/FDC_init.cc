// $Id$

#include "DEventLoop.h"
#include "DFactory_DFDCHit.h"
#include "DFactory_DFDCTruth.h"
#include "DFactory_DFDCPseudo.h"

derror_t FDC_init(DEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new DFactory_DFDCHit());
	loop->AddFactory(new DFactory_DFDCTruth());
	loop->AddFactory(new DFactory_DFDCPseudo());	

	return NOERROR;
}
