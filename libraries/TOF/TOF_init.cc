// $Id$

#include "DEvent.h"
#include "DFactory_DTOFHit.h"

derror_t TOF_init(DEvent *event)
{
	/// Create and register TOF data factories
	event->AddFactory(new DFactory_DTOFHit());

	return NOERROR;
}
