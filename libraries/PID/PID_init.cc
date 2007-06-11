// $Id: PID_init.cc 2433 2007-04-07 14:57:32Z kornicer $

#include "JANA/JEventLoop.h"
#include "DPhoton_factory.h"
#include "DPi0_factory.h"

#define UC_CLUSTERIZER

jerror_t PID_init(JEventLoop *loop)
{
	/// Create and register PID data factories
	loop->AddFactory(new DPhoton_factory());
	loop->AddFactory(new DPi0_factory());

	return NOERROR;
}
