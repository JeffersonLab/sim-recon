// $Id$

#include "JANA/JEventLoop.h"
#include "DTRIGGER_factory.h"

jerror_t TRIGGER_init(JEventLoop *loop)
{
	/// Create and register TRIGGER data factories
	loop->AddFactory(new DTRIGGER_factory());

	return NOERROR;
}
