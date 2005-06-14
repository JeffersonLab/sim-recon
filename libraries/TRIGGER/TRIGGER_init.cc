// $Id$

#include "DEventLoop.h"
#include "DFactory_DTRIGGER.h"

derror_t TRIGGER_init(DEventLoop *loop)
{
	/// Create and register TRIGGER data factories
	loop->AddFactory(new DFactory_DTRIGGER());

	return NOERROR;
}
