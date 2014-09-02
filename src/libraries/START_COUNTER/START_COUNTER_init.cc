// $Id: CDC_init.cc 2151 2006-10-23 18:01:37Z davidl $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DSCTruthHit.h"
#include "DSCDigiHit.h"
#include "DSCHit_factory.h"
#include "DSCTDCDigiHit.h"

jerror_t START_COUNTER_init(JEventLoop *loop)
{
	/// Create and register Start Counter data factories
	loop->AddFactory(new JFactory<DSCDigiHit>());
	loop->AddFactory(new JFactory<DSCTDCDigiHit>());
	loop->AddFactory(new DSCHit_factory());
	loop->AddFactory(new JFactory<DSCHit>("TRUTH"));
	loop->AddFactory(new JFactory<DSCTruthHit>());

	return NOERROR;
}
