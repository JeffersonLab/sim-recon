// $Id$

#include "DEvent.h"
#include "DFactory_TRACKINGHits.h"
#include "DFactory_MCCheatHits.h"
#include "DFactory_MCTrackCandidates.h"

derror_t TRACKING_init(DEvent *event)
{
	/// Create and register TRACKING data factories
	event->AddFactory(new DFactory_TRACKINGHits(event));
	event->AddFactory(new DFactory_MCCheatHits(event));
	event->AddFactory(new DFactory_MCTrackCandidates(event));

	return NOERROR;
}
