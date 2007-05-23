// $Id: PID_init.cc 2433 2007-04-07 14:57:32Z kornicer $

#include "JANA/JEventLoop.h"
#include "DPIDPhoton_factory.h"
#include "DPIDPi0_factory.h"

#define UC_CLUSTERIZER

jerror_t PID_init(JEventLoop *loop)
{
	/// Create and register PID data factories
	loop->AddFactory(new DPIDPhoton_factory());
	loop->AddFactory(new DPIDPi0_factory());

	return NOERROR;
}
