//***********************************************
// FDC_init.cc: defines a function called by
// the framework to register the FDC factories.
// Author: Craig Bookwalter
// Date: Apr 2006
//***********************************************

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DFDCHit_factory.h"
#include "DFDCPseudo_factory.h"
#include "DFDCCathodeCluster_factory.h"
#include "DFDCSegment_factory.h"
#include "DFDCIntersection_factory.h"
#include "DFDCPseudo_factory_WIRESONLY.h"

#include <FDC/DFDCCathodeDigiHit.h>
#include <FDC/DFDCWireDigiHit.h>

jerror_t FDC_init(JEventLoop *loop)
{
	/// Create and register FDC data factories
	loop->AddFactory(new JFactory<DFDCCathodeDigiHit>());
	loop->AddFactory(new JFactory<DFDCWireDigiHit>());
	loop->AddFactory(new DFDCHit_factory());
	loop->AddFactory(new JFactory<DFDCHit>("TRUTH"));
	loop->AddFactory(new JFactory<DFDCHit>("CALIB"));
	loop->AddFactory(new DFDCPseudo_factory());
	loop->AddFactory(new DFDCCathodeCluster_factory());
	loop->AddFactory(new DFDCSegment_factory());
	loop->AddFactory(new DFDCIntersection_factory());
	loop->AddFactory(new DFDCPseudo_factory_WIRESONLY());

	return NOERROR;
}
