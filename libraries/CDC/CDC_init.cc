// $Id$

#include "DEvent.h"
#include "DFactory_DCDCHit.h"

derror_t CDC_init(DEvent *event)
{
	/// Create and register CDC data factories
	event->AddFactory(new DFactory_DCDCHit());

	return NOERROR;
}
