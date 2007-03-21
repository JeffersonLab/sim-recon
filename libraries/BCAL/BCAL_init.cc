// $Id$

#include "JANA/JEventLoop.h"
#include "DBCALHit_factory.h"
#include "DHDDMBCALHit_factory.h"
#include "DBCALMCResponse_factory.h"
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory.h"
#include "DBCALShower_factory_SIMPLE.h"

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DHDDMBCALHit_factory());
	loop->AddFactory(new DBCALMCResponse_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory());
	loop->AddFactory(new DBCALShower_factory_SIMPLE());

	return NOERROR;
}
