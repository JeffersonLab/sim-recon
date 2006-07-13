// $Id$

#include "JANA/JEventLoop.h"
#include "DTAGGERHit_factory.h"

jerror_t TAGGER_init(JEventLoop *loop)
{
	/// Create and register TAGGER data factories
	loop->AddFactory(new DTAGGERHit_factory());

	return NOERROR;
}
