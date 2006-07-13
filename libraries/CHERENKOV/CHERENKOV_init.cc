// $Id$

#include "JANA/JEventLoop.h"
#include "DCHERENKOVHit_factory.h"

jerror_t CHERENKOV_init(JEventLoop *loop)
{
	/// Create and register CHERENKOV data factories
	loop->AddFactory(new DCHERENKOVHit_factory());

	return NOERROR;
}
