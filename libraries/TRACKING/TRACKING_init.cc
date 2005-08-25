// $Id$

#include "DEventLoop.h"
#include "DFactory_DTrack.h"
#include "DFactory_DTrackCandidate.h"
#include "DFactory_DTrackHit.h"
#include "DFactory_DTrackHit_MC.h"
#include "DFactory_DMCTrackHit.h"
#include "DFactory_DTrackEfficiency.h"
#include "DFactory_DMCThrown.h"

derror_t TRACKING_init(DEventLoop *loop)
{
	/// Create and register TRACKING data factories
	loop->AddFactory(new DFactory_DTrack());
	loop->AddFactory(new DFactory_DTrackCandidate());
	loop->AddFactory(new DFactory_DTrackHit());
	loop->AddFactory(new DFactory_DTrackHit_MC());
	loop->AddFactory(new DFactory_DMCTrackHit());
	loop->AddFactory(new DFactory_DTrackEfficiency());
	loop->AddFactory(new DFactory_DMCThrown());

	return NOERROR;
}
