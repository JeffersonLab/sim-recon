// $Id$

#include "DEvent.h"
#include "DFactory_DFCALHit.h"
#include "DFactory_DFCALShowers.h"

derror_t FCAL_init(DEvent *event)
{
	/// Create and register FCAL data factories
	event->AddFactory(new DFactory_DFCALHit());
	event->AddFactory(new DFactory_DFCALShowers());

	return NOERROR;
}
