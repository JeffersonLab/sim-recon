// $Id$

#include "DEvent.h"
#include "DFactory_TRACKINGHits.h"

derror_t TRACKING_init(DEvent *event)
{
	/// Create and register TRACKING data factories
	event->AddFactory(new DFactory_TRACKINGHits(event));

	return NOERROR;
}
