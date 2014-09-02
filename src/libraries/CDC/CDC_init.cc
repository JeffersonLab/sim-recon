// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DCDCDigiHit.h"
#include "DCDCHit_factory.h"
#include "DCDCTrackHit_factory.h"


jerror_t CDC_init(JEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new JFactory<DCDCDigiHit>());
	loop->AddFactory(new DCDCHit_factory());
	loop->AddFactory(new JFactory<DCDCHit>("TRUTH"));
	loop->AddFactory(new DCDCTrackHit_factory());

	return NOERROR;
}
