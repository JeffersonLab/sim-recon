// $Id$

#include "DEventLoop.h"
#include "DFactory_DFCALHit.h"
#include "DFactory_DFCALShower.h"

derror_t FCAL_init(DEventLoop *loop)
{
	/// Create and register FCAL data factories
	loop->AddFactory(new DFactory_DFCALHit());
	loop->AddFactory(new DFactory_DFCALShower());

	return NOERROR;
}
