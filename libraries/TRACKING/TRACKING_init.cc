// $Id$

#include "JANA/JEventLoop.h"
#include "DTrack_factory.h"
#include "DTrackCandidate_factory.h"
#include "DTrackCandidate_factory_THROWN.h"
#include "DTrackCandidate_factory_CDC.h"
#include "DTrackCandidate_factory_FDC.h"
#include "DTrackCandidate_factory_FDCCathodes.h"
#include "DTrackCandidate_factory_FDCpseudo.h"
#include "DTrackHit_factory.h"
#include "DTrackHit_factory_MC.h"
#include "DMCTrackHit_factory.h"
#include "DMCThrown_factory.h"
#include "DMCTrajectoryPoint_factory.h"
#include "DTrack_factory_THROWN.h"
#include "DTrack_factory_ALT1.h"
#include "DTrack_factory_ALT2.h"


jerror_t TRACKING_init(JEventLoop *loop)
{
	/// Create and register TRACKING data factories
	loop->AddFactory(new DTrack_factory_ALT1());
	loop->AddFactory(new DTrack_factory_ALT2());
	loop->AddFactory(new DTrack_factory());
	loop->AddFactory(new DTrackCandidate_factory());
	loop->AddFactory(new DTrackCandidate_factory_CDC());
	loop->AddFactory(new DTrackCandidate_factory_FDC());
	loop->AddFactory(new DTrackCandidate_factory_FDCCathodes());
	loop->AddFactory(new DTrackCandidate_factory_FDCpseudo());
	loop->AddFactory(new DTrackCandidate_factory_THROWN());
	loop->AddFactory(new DTrackHit_factory());
	loop->AddFactory(new DTrackHit_factory_MC());
	loop->AddFactory(new DMCTrackHit_factory());
	loop->AddFactory(new DMCThrown_factory());
	loop->AddFactory(new DMCTrajectoryPoint_factory());
	loop->AddFactory(new DTrack_factory_THROWN());

	return NOERROR;
}
