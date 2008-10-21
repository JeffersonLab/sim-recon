// $Id$

#include <JANA/JEventLoop.h>
#include "DBCALMCResponse_factory.h"
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory.h"
#include "DBCALPhoton_factory.h"

#include "DBCALHit.h"
#include "DBCALTruthShower.h"
#include "DHDDMBCALHit.h"
typedef JFactory<DBCALHit> DBCALHit_factory;
typedef JFactory<DBCALTruthShower> DBCALTruthShower_factory;
typedef JFactory<DHDDMBCALHit> DMCBCALHit_factory;

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DBCALMCResponse_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory());
	loop->AddFactory(new DBCALTruthShower_factory());
	loop->AddFactory(new DBCALPhoton_factory());
	loop->AddFactory(new DMCBCALHit_factory());
    
	return NOERROR;
}
