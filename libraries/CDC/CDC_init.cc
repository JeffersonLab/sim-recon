// $Id$

#include "DEvent.h"
#include "DFactory_CDCHits.h"
#include "DFactory_CDCClusters.h"

derror_t CDC_init(DEvent *event)
{
	/// Create and register CDC data factories
	event->AddFactory(new DFactory_CDCHits(event));
	event->AddFactory(new DFactory_CDCClusters(event));

	return NOERROR;
}
