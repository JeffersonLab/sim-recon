//***********************************************
// FDC_init.cc: defines a function called by
// the framework to register the FDC factories.
// Author: Craig Bookwalter
// Date: Apr 2006
//***********************************************

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DFDCPseudo_factory.h"
#include "DFDCCathodeCluster_factory.h"
#include "DFDCSegment_factory.h"
#include "DFDCIntersection_factory.h"
#include "DFDCPseudo_factory_WIRESONLY.h"

#include "DFDCHit.h"
typedef JFactory<DFDCHit> DFDCHit_factory;

jerror_t FDC_init(JEventLoop *loop)
{
	jout<<"Registering FDC factories"<<endl;
  
	/// Create and register FDC data factories
	loop->AddFactory(new DFDCHit_factory());
	loop->AddFactory(new DFDCHit_factory("TRUTH"));
	loop->AddFactory(new DFDCPseudo_factory());
	loop->AddFactory(new DFDCCathodeCluster_factory());
	loop->AddFactory(new DFDCSegment_factory());
	loop->AddFactory(new DFDCIntersection_factory());
	loop->AddFactory(new DFDCPseudo_factory_WIRESONLY());

	return NOERROR;
}
