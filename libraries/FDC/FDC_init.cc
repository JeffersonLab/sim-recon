//***********************************************
// FDC_init.cc: defines a function called by
// the framework to register the FDC factories.
// Author: Craig Bookwalter
// Date: Apr 2006
//***********************************************

#include "JANA/JEventLoop.h"
#include "DFDCHit_factory.h"
#include "DFDCTruth_factory.h"
#include "DFDCPseudo_factory.h"
#include "DFDCCathodeCluster_factory.h"
#include "DFDCSegment_factory.h"

jerror_t FDC_init(JEventLoop *loop)
{
  printf("Registering FDC factories");
  
	/// Create and register FDC data factories
	loop->AddFactory(new DFDCHit_factory());
	loop->AddFactory(new DFDCTruth_factory());
	loop->AddFactory(new DFDCPseudo_factory());
	loop->AddFactory(new DFDCCathodeCluster_factory());
	loop->AddFactory(new DFDCSegment_factory());

	return NOERROR;
}
