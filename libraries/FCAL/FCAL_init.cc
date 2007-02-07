// $Id$

#include "JANA/JEventLoop.h"
#include "DFCALHit_factory.h"
#include "DFCALShower_factory.h"
#include "DFCALCluster_factory.h"
#include "DFCALPhoton_factory.h"
#include "DFCALGeometry_factory.h"
#include "DFCALTruthShower_factory.h"
#include "DFCALMCResponse_factory.h"

#define UC_CLUSTERIZER

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFCALHit_factory());
#ifdef UC_CLUSTERIZER
	loop->AddFactory(new DFCALCluster_factory());
	loop->AddFactory(new DFCALPhoton_factory());
#else
	loop->AddFactory(new DFCALShower_factory());
#endif
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new DFCALTruthShower_factory());
	loop->AddFactory(new DFCALMCResponse_factory());

	return NOERROR;
}
