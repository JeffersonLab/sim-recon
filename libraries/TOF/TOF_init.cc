// $Id$

#include "JANA/JEventLoop.h"
#include "DTOFHit_factory.h"
#include "DTOFGeometry_factory.h"
#include "DTOFMCResponse_factory.h"
#include "DTOFTruth_factory.h"
#include "DHDDMTOFHit_factory.h"
#include "DTOFPoint_factory.h"

jerror_t TOF_init(JEventLoop *loop)
{
	/// Create and register TOF data factories
	loop->AddFactory(new DTOFMCResponse_factory());
	loop->AddFactory(new DTOFHit_factory());
	loop->AddFactory(new DTOFGeometry_factory());
	loop->AddFactory(new DTOFTruth_factory());
	loop->AddFactory(new DHDDMTOFHit_factory());
	loop->AddFactory(new DTOFPoint_factory());

	return NOERROR;
}
