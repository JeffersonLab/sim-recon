// $Id$

#include "DEvent.h"
#include "DFactory_TRIGGERHits.h"

derror_t TRIGGER_init(DEvent *event)
{
	/// Create and register TRIGGER data factories
	event->AddFactory(new DFactory_TRIGGERHits(event));

	return NOERROR;
}
