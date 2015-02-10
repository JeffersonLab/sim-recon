// $Id: CDC_init.cc 2151 2006-10-23 18:01:37Z davidl $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTPOLTruthHit.h"
#include "DTPOLSectorDigiHit.h"
#include "DTPOLRingDigiHit.h"
#include "DTPOLHit_factory.h"

jerror_t TPOL_init(JEventLoop *loop)
{
	/// Create and register TPOL data factories
	loop->AddFactory(new JFactory<DTPOLSectorDigiHit>());
	loop->AddFactory(new JFactory<DTPOLRingDigiHit>());
	loop->AddFactory(new DTPOLHit_factory());
	loop->AddFactory(new JFactory<DTPOLHit>("TRUTH"));
	loop->AddFactory(new JFactory<DTPOLTruthHit>());

	return NOERROR;
}
