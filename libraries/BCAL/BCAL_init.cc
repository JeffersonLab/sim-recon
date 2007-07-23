// $Id$

#include "JANA/JEventLoop.h"
#include "BCAL/DBCALHit_factory.h"
#include "BCAL/DBCALTruthShower_factory.h"
#include "BCAL/DBCALMCResponse_factory.h"
#include "BCAL/DBCALGeometry_factory.h"
#include "BCAL/DBCALShower_factory.h"
#include "BCAL/DBCALShower_factory_SIMPLE.h"
#include "BCAL/DBCALShower_factory_IU.h"
#include "BCAL/DBCALPhoton_factory.h"

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DBCALMCResponse_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory());
	loop->AddFactory(new DBCALShower_factory_SIMPLE());
    loop->AddFactory(new DBCALShower_factory_IU());
	loop->AddFactory(new DBCALTruthShower_factory());
    loop->AddFactory(new DBCALPhoton_factory());
    
	return NOERROR;
}
