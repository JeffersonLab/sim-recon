// $Id$

#include "DEvent.h"
#include "DFactory_FDCHits.h"
#include "DFactory_FDCClusters.h"

derror_t FDC_init(DEvent *event)
{
	/// Create and register FDC data factories
	event->AddFactory(new DFactory_FDCHits(event));
	event->AddFactory(new DFactory_FDCClusters(event));

	return NOERROR;
}
