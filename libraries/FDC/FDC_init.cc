// $Id$

#include "DEvent.h"
#include "DFactory_DFDCHit.h"

derror_t FDC_init(DEvent *event)
{
	/// Create and register FDC data factories
	event->AddFactory(new DFactory_DFDCHit());

	return NOERROR;
}
