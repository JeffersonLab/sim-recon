// $Id$

#include "JANA/JEventLoop.h"
#include "DCDCHit_factory.h"

jerror_t CDC_init(JEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new DCDCHit_factory());

	return NOERROR;
}
