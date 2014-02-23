// $Id: CDC_init.cc 2151 2006-10-23 18:01:37Z davidl $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DSCTruthHit.h"
#include "DSCDigiHit.h"
#include "DSCHit_factory.h"
#include "DSCTDCDigiHit.h"
typedef JFactory<DSCDigiHit> DSCDigiHit_factory;
typedef JFactory<DSCTDCDigiHit> DSCTDCDigiHit_factory;
typedef JFactory<DSCTruthHit> DSCTruthHit_factory;

jerror_t START_COUNTER_init(JEventLoop *loop)
{
	/// Create and register Start Counter data factories
	loop->AddFactory(new DSCDigiHit_factory());
	loop->AddFactory(new DSCTDCDigiHit_factory());
	loop->AddFactory(new DSCHit_factory());

	loop->AddFactory(new DSCTruthHit_factory());

	return NOERROR;
}
