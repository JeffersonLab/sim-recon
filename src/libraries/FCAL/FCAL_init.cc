// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DFCALCluster_factory.h"
#include "DFCALGeometry_factory.h"
#include "DFCALShower_factory.h"
#include "DFCALTruthShower.h"
#include "DFCALHit.h"
typedef JFactory<DFCALTruthShower> DFCALTruthShower_factory;
typedef JFactory<DFCALHit> DFCALHit_factory;

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFCALHit_factory());
	loop->AddFactory(new DFCALHit_factory("TRUTH"));
	loop->AddFactory(new DFCALCluster_factory());
	loop->AddFactory(new DFCALShower_factory());
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new DFCALTruthShower_factory());

	return NOERROR;
}
