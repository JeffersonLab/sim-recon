// $Id: PID_init.cc 2433 2007-04-07 14:57:32Z kornicer $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DBeamPhoton_factory.h"
#include "DBeamPhoton_factory_TRUTH.h"
#include "DBeamPhoton_factory_TAGGEDMCGEN.h"
#include "DBeamPhoton_factory_MCGEN.h"
#include "DParticleID_factory.h"
#include "DParticleID_factory_PID1.h"
#include "DChargedTrack_factory.h"
#include "DChargedTrack_factory_PreSelect.h"
#include "DChargedTrackHypothesis_factory.h"
#include "DNeutralParticle_factory.h"
#include "DNeutralParticle_factory_PreSelect.h"
#include "DNeutralParticleHypothesis_factory.h"
#include "DNeutralShower_factory.h"
#include "DNeutralShower_factory_PreSelect.h"
#include "DVertex_factory.h"
#include "DVertex_factory_THROWN.h"
#include "DEventRFBunch_factory.h"
#include "DEventRFBunch_factory_Thrown.h"
#include "DEventRFBunch_factory_Calibrations.h"
#include "DDetectorMatches_factory.h"
#include "DDetectorMatches_factory_WireBased.h"
#include "DMCThrown_factory_FinalState.h"
#include "DMCThrown_factory_Decaying.h"
#include "DMCThrown_factory_Primary.h"

#include "DMCReaction.h"

#define UC_CLUSTERIZER

jerror_t PID_init(JEventLoop *loop)
{
	/// Create and register PID data factories
	loop->AddFactory(new JFactory<DMCReaction>());
	loop->AddFactory(new DBeamPhoton_factory);
	loop->AddFactory(new DBeamPhoton_factory_TRUTH);
	loop->AddFactory(new DBeamPhoton_factory_TAGGEDMCGEN);
	loop->AddFactory(new DBeamPhoton_factory_MCGEN);
	loop->AddFactory(new DParticleID_factory);
	loop->AddFactory(new DParticleID_factory_PID1);
	loop->AddFactory(new DChargedTrack_factory);
	loop->AddFactory(new DChargedTrack_factory_PreSelect);
	loop->AddFactory(new DChargedTrackHypothesis_factory);
	loop->AddFactory(new DNeutralParticle_factory);
	loop->AddFactory(new DNeutralParticle_factory_PreSelect);
	loop->AddFactory(new DNeutralParticleHypothesis_factory);
	loop->AddFactory(new DNeutralShower_factory);
	loop->AddFactory(new DNeutralShower_factory_PreSelect);
	loop->AddFactory(new DVertex_factory);
	loop->AddFactory(new DVertex_factory_THROWN);
	loop->AddFactory(new DEventRFBunch_factory);
	loop->AddFactory(new DEventRFBunch_factory_Thrown);
	loop->AddFactory(new DEventRFBunch_factory_Calibrations);
	loop->AddFactory(new DDetectorMatches_factory);
	loop->AddFactory(new DDetectorMatches_factory_WireBased);
	loop->AddFactory(new DMCThrown_factory_FinalState);
	loop->AddFactory(new DMCThrown_factory_Decaying);
	loop->AddFactory(new DMCThrown_factory_Primary);

	return NOERROR;
}
