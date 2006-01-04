// $Id$

#include "DEventLoop.h"
#include "DFactory_DFCALHit.h"
#include "DFactory_DFCALShower.h"
#include "DFactory_DFCALGeometry.h"
#include "DFactory_DFCALTruthShower.h"
#include "DFactory_DFCALMCResponse.h"

derror_t FCAL_init(DEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFactory_DFCALHit());
	loop->AddFactory(new DFactory_DFCALShower());
	loop->AddFactory(new DFactory_DFCALGeometry());
	loop->AddFactory(new DFactory_DFCALTruthShower());
	loop->AddFactory(new DFactory_DFCALMCResponse());

	return NOERROR;
}
