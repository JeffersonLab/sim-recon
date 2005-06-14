// $Id$

#include "DEventLoop.h"
#include "DFactory_DTAGGERHit.h"

derror_t TAGGER_init(DEventLoop *loop)
{
	/// Create and register TAGGER data factories
	loop->AddFactory(new DFactory_DTAGGERHit());

	return NOERROR;
}
