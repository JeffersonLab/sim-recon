// $Id$

#include "DEvent.h"
#include "DFactory_DUPVHit.h"

derror_t UPV_init(DEvent *event)
{
	/// Create and register UPV data factories
	event->AddFactory(new DFactory_DUPVHit());

	return NOERROR;
}
