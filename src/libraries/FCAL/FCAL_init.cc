// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DFCALCluster_factory.h"
#include "DFCALGeometry_factory.h"
#include "DFCALShower_factory.h"
#include "DFCALTruthShower.h"
#include "DFCALDigiHit.h"
#include "DFCALHit_factory.h"

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new JFactory<DFCALDigiHit>());
	loop->AddFactory(new DFCALHit_factory());
	loop->AddFactory(new JFactory<DFCALHit>("TRUTH"));
	loop->AddFactory(new DFCALCluster_factory());
	loop->AddFactory(new DFCALShower_factory());
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new JFactory<DFCALTruthShower>());

	return NOERROR;
}
