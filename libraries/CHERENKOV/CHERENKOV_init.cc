// $Id$

#include "DEvent.h"
#include "DFactory_CHERENKOVHits.h"

derror_t CHERENKOV_init(DEvent *event)
{
	/// Create and register CHERENKOV data factories
	event->AddFactory(new DFactory_CHERENKOVHits(event));

	return NOERROR;
}
