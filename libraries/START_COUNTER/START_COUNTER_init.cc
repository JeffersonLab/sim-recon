// $Id: CDC_init.cc 2151 2006-10-23 18:01:37Z davidl $

#include "JANA/JEventLoop.h"
#include "DSCTruthHit_factory.h"
#include "DSCHit_factory.h"

jerror_t START_COUNTER_init(JEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new DSCTruthHit_factory());
	loop->AddFactory(new DSCHit_factory());

	return NOERROR;
}
