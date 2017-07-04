#include "PID/DKinematicData.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"

//Declare thread_local resource pools
thread_local DResourcePool<DKinematicData::DKinematicInfo> DKinematicData::dResourcePool_KinematicInfo;

