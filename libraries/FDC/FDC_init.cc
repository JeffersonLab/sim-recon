//***********************************************
// FDC_init.cc: defines a function called by
// the framework to register the FDC factories.
// Author: Craig Bookwalter
// Date: Apr 2006
//***********************************************

#include "DEventLoop.h"
#include "DFactory_DFDCHit.h"
#include "DFactory_DFDCTruth.h"
#include "DFactory_DFDCPseudo.h"
#include "DFactory_DFDCCathodeCluster.h"

derror_t FDC_init(DEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new DFactory_DFDCHit());
	loop->AddFactory(new DFactory_DFDCTruth());
	loop->AddFactory(new DFactory_DFDCPseudo());
	loop->AddFactory(new DFactory_DFDCCathodeCluster());

	return NOERROR;
}
