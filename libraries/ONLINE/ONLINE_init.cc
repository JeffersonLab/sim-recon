// $Id: ONLINE_init.cc 2167 2006-11-06 15:53:12Z davidl $

#include "JANA/JEventLoop.h"
#include "DFADC_factory.h"

jerror_t ONLINE_init(JEventLoop *loop)
{
	/// Create and register ONLINE data factories
	loop->AddFactory(new DFADC_factory());

	return NOERROR;
}
