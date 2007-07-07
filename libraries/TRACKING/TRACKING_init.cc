// $Id$

#include "JANA/JEventLoop.h"
#include "DTrack_factory.h"
#include "DTrackCandidate_factory.h"
#include "DTrackCandidate_factory_THROWN.h"
#include "DTrackHit_factory.h"
#include "DTrackHit_factory_MC.h"
#include "DMCTrackHit_factory.h"
#include "DTrackEfficiency_factory.h"
#include "DMCThrown_factory.h"
#include "DMCTrajectoryPoint_factory.h"
#include "DTrackLinker_factory.h"


jerror_t TRACKING_init(JEventLoop *loop)
{
	/// Create and register TRACKING data factories
	loop->AddFactory(new DTrack_factory());
	loop->AddFactory(new DTrackCandidate_factory());
	loop->AddFactory(new DTrackCandidate_factory_THROWN());
	loop->AddFactory(new DTrackHit_factory());
	loop->AddFactory(new DTrackHit_factory_MC());
	loop->AddFactory(new DMCTrackHit_factory());
	loop->AddFactory(new DTrackEfficiency_factory());
	loop->AddFactory(new DMCThrown_factory());
	loop->AddFactory(new DMCTrajectoryPoint_factory());
	loop->AddFactory(new DTrackLinker_factory());

	return NOERROR;
}
