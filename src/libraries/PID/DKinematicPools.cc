#include "PID/DKinematicData.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"

//Declare thread_local resource pools
thread_local DResourcePool<DKinematicData::DKinematicInfo> DKinematicData::dResourcePool_KinematicInfo;
thread_local DResourcePool<DChargedTrackHypothesis::DTimingInfo> DChargedTrackHypothesis::dResourcePool_TimingInfo;
thread_local DResourcePool<DChargedTrackHypothesis::DTrackingInfo> DChargedTrackHypothesis::dResourcePool_TrackingInfo;
thread_local DResourcePool<DNeutralParticleHypothesis::DTimingInfo> DNeutralParticleHypothesis::dResourcePool_TimingInfo;

