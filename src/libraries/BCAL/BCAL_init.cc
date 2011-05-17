// $Id$

#include <JANA/JEventLoop.h>
#include "DBCALGeometry_factory.h"
#include "DBCALShower_factory_KLOE.h"
#include "DBCALShower_factory.h"
#include "DBCALCluster_factory.h"
#include "DBCALPhoton_factory.h"
#include "DBCALHit.h"
#include "DBCALSiPMHit.h"
#include "DBCALTruthCell.h"

#include "DBCALTruthShower.h"

// These come from the event source, not from any algorithm
typedef JFactory<DBCALHit> DBCALHit_factory;
typedef JFactory<DBCALSiPMHit> DBCALSiPMHit_factory;
typedef JFactory<DBCALTruthShower> DBCALTruthShower_factory;
typedef JFactory<DBCALTruthCell> DBCALTruthCell_factory;

jerror_t BCAL_init(JEventLoop *loop)
{
	/// Create and register BCAL data factories
	loop->AddFactory(new DBCALHit_factory());
	loop->AddFactory(new DBCALSiPMHit_factory());
	loop->AddFactory(new DBCALGeometry_factory());
	loop->AddFactory(new DBCALShower_factory_KLOE());
  loop->AddFactory(new DBCALShower_factory());
  loop->AddFactory(new DBCALCluster_factory());
	loop->AddFactory(new DBCALTruthShower_factory());
	loop->AddFactory(new DBCALPhoton_factory());
	loop->AddFactory(new DBCALTruthCell_factory());
    
	return NOERROR;
}
