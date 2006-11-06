// $Id$

#include "JANA/JEventLoop.h"
#include "DUPVHit_factory.h"
#include "DUPVTruthHit_factory.h"

jerror_t UPV_init(JEventLoop *loop)
{
	/// Create and register UPV data factories
	loop->AddFactory(new DUPVHit_factory());
	loop->AddFactory(new DUPVTruthHit_factory());

	return NOERROR;
}
