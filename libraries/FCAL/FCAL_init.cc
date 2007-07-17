// $Id$

#include "JANA/JEventLoop.h"
#include "FCAL/DFCALHit_factory.h"
#include "FCAL/DFCALCluster_factory.h"
#include "FCAL/DFCALPhoton_factory.h"
#include "FCAL/DFCALGeometry_factory.h"
#include "FCAL/DFCALTruthShower_factory.h"
#include "FCAL/DFCALMCResponse_factory.h"

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFCALHit_factory());
	loop->AddFactory(new DFCALCluster_factory());
	loop->AddFactory(new DFCALPhoton_factory());
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new DFCALTruthShower_factory());
	loop->AddFactory(new DFCALMCResponse_factory());

	return NOERROR;
}
