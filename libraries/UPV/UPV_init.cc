// $Id$

#include "DEvent.h"
#include "DFactory_UPVHits.h"

derror_t UPV_init(DEvent *event)
{
	/// Create and register UPV data factories
	event->AddFactory(new DFactory_UPVHits(event));

	return NOERROR;
}
