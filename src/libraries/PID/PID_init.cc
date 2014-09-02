// $Id: PID_init.cc 2433 2007-04-07 14:57:32Z kornicer $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DBeamPhoton_factory.h"
#include "DBeamPhoton_factory_TRUTH.h"
#include "DParticleID_factory.h"
#include "DParticleID_factory_PID1.h"
#include "DChargedTrack_factory.h"
#include "DChargedTrackHypothesis_factory.h"
#include "DNeutralParticle_factory.h"
#include "DNeutralParticleHypothesis_factory.h"
#include "DNeutralShower_factory.h"
#include "DVertex_factory.h"
#include "DVertex_factory_THROWN.h"
#include "DEventRFBunch_factory.h"
#include "DEventRFBunch_factory_Thrown.h"
#include "DDetectorMatches_factory.h"
#include "DMCThrown_factory_FinalState.h"
#include "DMCThrown_factory_Decaying.h"
#include "DMCThrown_factory_Primary.h"

#include "DRFTime.h"
#include "DBeamPhoton.h"
#include "DMCReaction.h"

#define UC_CLUSTERIZER

jerror_t PID_init(JEventLoop *loop)
{
	/// Create and register PID data factories
	loop->AddFactory(new JFactory<DMCReaction>());
	loop->AddFactory(new JFactory<DRFTime>());
	loop->AddFactory(new JFactory<DBeamPhoton>("MCGEN"));
	loop->AddFactory(new DBeamPhoton_factory);
	loop->AddFactory(new DBeamPhoton_factory_TRUTH);
	loop->AddFactory(new DParticleID_factory);
	loop->AddFactory(new DParticleID_factory_PID1);
	loop->AddFactory(new DChargedTrack_factory);
	loop->AddFactory(new DChargedTrackHypothesis_factory);
	loop->AddFactory(new DNeutralParticle_factory);
	loop->AddFactory(new DNeutralParticleHypothesis_factory);
	loop->AddFactory(new DNeutralShower_factory);
	loop->AddFactory(new DVertex_factory);
	loop->AddFactory(new DVertex_factory_THROWN);
	loop->AddFactory(new DEventRFBunch_factory);
	loop->AddFactory(new DEventRFBunch_factory_Thrown);
	loop->AddFactory(new DDetectorMatches_factory);
	loop->AddFactory(new DMCThrown_factory_FinalState);
	loop->AddFactory(new DMCThrown_factory_Decaying);
	loop->AddFactory(new DMCThrown_factory_Primary);

	return NOERROR;
}
