// $Id$

#include "DEvent.h"
#include "DFactory_BCALHits.h"

derror_t BCAL_init(DEvent *event)
{
	/// Create and register BCAL data factories
	event->AddFactory(new DFactory_BCALHits(event));

	return NOERROR;
}
