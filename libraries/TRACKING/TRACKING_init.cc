// $Id$

#include "DEvent.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCTrackCandidate.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DMCReconstructed.h"

derror_t TRACKING_init(DEvent *event)
{
	/// Create and register TRACKING data factories
	event->AddFactory(new DFactory_DMCCheatHit());
	event->AddFactory(new DFactory_DMCTrackCandidate());
	event->AddFactory(new DFactory_DMCThrown());
	event->AddFactory(new DFactory_DMCReconstructed());

	return NOERROR;
}
