// $Id$

#include "JANA/JEventLoop.h"
#include "DFCALHit_factory.h"
#include "DFCALShower_factory.h"
#include "DFCALGeometry_factory.h"
#include "DFCALTruthShower_factory.h"
#include "DFCALMCResponse_factory.h"

jerror_t FCAL_init(JEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFCALHit_factory());
	loop->AddFactory(new DFCALShower_factory());
	loop->AddFactory(new DFCALGeometry_factory());
	loop->AddFactory(new DFCALTruthShower_factory());
	loop->AddFactory(new DFCALMCResponse_factory());

	return NOERROR;
}
