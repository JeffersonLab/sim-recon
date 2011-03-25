// $Id$

#include <JANA/JEventLoop.h>
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory_KLOE.h"
#include "DBCALShower_factory.h"
#include "DBCALCluster_factory.h"
#include "DBCALPhoton_factory.h"
#include "DBCALHit.h"

#include "DBCALTruthShower.h"

typedef JFactory<DBCALHit> DBCALHit_factory;
typedef JFactory<DBCALTruthShower> DBCALTruthShower_factory;

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory_KLOE());
  loop->AddFactory(new DBCALShower_factory());
  loop->AddFactory(new DBCALCluster_factory());
	loop->AddFactory(new DBCALTruthShower_factory());
	loop->AddFactory(new DBCALPhoton_factory());
    
	return NOERROR;
}
