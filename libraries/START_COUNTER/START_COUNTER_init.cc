// $Id: CDC_init.cc 2151 2006-10-23 18:01:37Z davidl $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DSCTruthHit.h"
#include "DSCHit.h"
typedef JFactory<DSCHit> DSCHit_factory;
typedef JFactory<DSCTruthHit> DSCTruthHit_factory;

jerror_t START_COUNTER_init(JEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new DSCTruthHit_factory());
	loop->AddFactory(new DSCHit_factory());

	return NOERROR;
}
