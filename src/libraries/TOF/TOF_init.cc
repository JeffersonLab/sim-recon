// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTOFHit_factory.h"
#include "DTOFMCResponse_factory.h"
#include "DTOFPoint_factory.h"
#include "DTOFHit_factory_MC.h"
#include "DTOFGeometry_factory.h"

#include "DHDDMTOFHit.h"
#include "DTOFTruth.h"
typedef JFactory<DHDDMTOFHit> DHDDMTOFHit_factory;
typedef JFactory<DTOFTruth> DTOFTruth_factory;

jerror_t TOF_init(JEventLoop *loop)
{
	/// Create and register TOF data factories
	loop->AddFactory(new DTOFMCResponse_factory());
	loop->AddFactory(new DTOFHit_factory());
	loop->AddFactory(new DTOFGeometry_factory());
	loop->AddFactory(new DTOFTruth_factory());
	loop->AddFactory(new DHDDMTOFHit_factory());
	loop->AddFactory(new DTOFPoint_factory());
	loop->AddFactory(new DTOFHit_factory_MC());

	return NOERROR;
}
