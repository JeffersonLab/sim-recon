// $Id$

#include "DEvent.h"
#include "DFactory_DFCALHit.h"

derror_t FCAL_init(DEvent *event)
{
	/// Create and register FCAL data factories
	event->AddFactory(new DFactory_DFCALHit());

	return NOERROR;
}
