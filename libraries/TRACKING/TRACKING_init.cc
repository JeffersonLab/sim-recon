// $Id$

#include "DEventLoop.h"
#include "DFactory_DMCCheatHit.h"
#include "DFactory_DMCTrackCandidate.h"
#include "DFactory_DMCThrown.h"
#include "DFactory_DMCReconstructed.h"
#include "DFactory_DMCTrackEfficiency.h"
#include "DFactory_DMCTrackCandidate_B.h"

derror_t TRACKING_init(DEventLoop *loop)
{
	/// Create and register TRACKING data factories
	loop->AddFactory(new DFactory_DMCCheatHit());
	loop->AddFactory(new DFactory_DMCTrackCandidate());
	loop->AddFactory(new DFactory_DMCThrown());
	loop->AddFactory(new DFactory_DMCReconstructed());
	loop->AddFactory(new DFactory_DMCTrackEfficiency());
	loop->AddFactory(new DFactory_DMCTrackCandidate_B());

	return NOERROR;
}
