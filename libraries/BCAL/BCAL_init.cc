// $Id$

#include "DEvent.h"
#include "DFactory_DBCALHit.h"

derror_t BCAL_init(DEvent *event)
{
	/// Create and register BCAL data factories
	event->AddFactory(new DFactory_DBCALHit());

	return NOERROR;
}
