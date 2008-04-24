// $Id$

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DFCALCluster_factory.h"
#include "DFCALGeometry_factory.h"
#include "DFCALHit_factory.h"
#include "DFCALMCResponse_factory.h"
#include "DFCALPhoton_factory.h"

#include "DFCALTruthShower.h"
#include "DMCFCALHit.h"
typedef JFactory<DFCALTruthShower> DFCALTruthShower_factory;
typedef JFactory<DMCFCALHit> DMCFCALHit_factory;

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFCALHit_factory());
	loop->AddFactory(new DFCALCluster_factory());
	loop->AddFactory(new DFCALPhoton_factory());
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new DFCALTruthShower_factory());
	loop->AddFactory(new DFCALMCResponse_factory());

	// This is just used as a container
	loop->AddFactory(new JFactory<DMCFCALHit>());

	return NOERROR;
}
