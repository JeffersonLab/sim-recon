// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DCDCHit.h"
#include "DCDCTrackHit_factory.h"

typedef JFactory<DCDCHit> DCDCHit_factory;

jerror_t CDC_init(JEventLoop *loop)
{
	/// Create and register CDC data factories
	loop->AddFactory(new DCDCHit_factory());
	loop->AddFactory(new DCDCHit_factory("TRUTH"));
	loop->AddFactory(new DCDCTrackHit_factory());

	return NOERROR;
}
