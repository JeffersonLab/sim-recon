// $Id: FCAL_init.cc 5342 2009-07-09 15:46:12Z davidl $

#include <JANA/JEventLoop.h>
using namespace jana;

#include <CCAL/DCCALTruthShower_factory.h>
#include <CCAL/DCCALGeometry_factory.h>

#include "DCCALHit.h"
typedef JFactory<DCCALHit> DCCALHit_factory;

jerror_t CCAL_init(JEventLoop *loop)
{
	/// Create and register CCAL data factories
	loop->AddFactory(new DCCALHit_factory());
	loop->AddFactory(new DCCALHit_factory("TRUTH"));
	loop->AddFactory(new DCCALTruthShower_factory());
	loop->AddFactory(new DCCALGeometry_factory());

	return NOERROR;
}
