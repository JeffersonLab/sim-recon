// $Id: PID_init.cc 2433 2007-04-07 14:57:32Z kornicer $

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DPhoton_factory.h"
#include "DPhoton_factory_THROWN.h"
#include "DTwoGammaFit_factory.h"
#include "DTwoGammaFit_factory_PI0.h"
#include "DTwoGammaFit_factory_ETA.h"
#include "DParticle_factory.h"
#include "DParticle_factory_THROWN.h"
#include "DParticle_factory_Kalman.h"
#include "DChargedTrack_factory.h"
#include "DChargedTrack_factory_Kalman.h"
#include "DChargedTruthMatch_factory.h"
#include "DVertex_factory.h"
#include "DPhysicsEvent_factory.h"

#include "DBeamPhoton.h"
typedef JFactory<DBeamPhoton> DBeamPhoton_factory;

#define UC_CLUSTERIZER

jerror_t PID_init(JEventLoop *loop)
{
	/// Create and register PID data factories
	loop->AddFactory(new DPhoton_factory());
//	loop->AddFactory(new DPi0_factory());
//	loop->AddFactory(new DTwoGammaFit_factory(0.135));
	loop->AddFactory(new DTwoGammaFit_factory_PI0);
	loop->AddFactory(new DTwoGammaFit_factory_ETA);
	loop->AddFactory(new DBeamPhoton_factory);
	loop->AddFactory(new DParticle_factory);
	loop->AddFactory(new DParticle_factory_THROWN);
	loop->AddFactory(new DParticle_factory_Kalman);
	loop->AddFactory(new DChargedTrack_factory);
	loop->AddFactory(new DChargedTrack_factory_Kalman);
	loop->AddFactory(new DChargedTruthMatch_factory);
	loop->AddFactory(new DPhoton_factory_THROWN);
	loop->AddFactory(new DVertex_factory);
	loop->AddFactory(new DPhysicsEvent_factory);

	return NOERROR;
}
