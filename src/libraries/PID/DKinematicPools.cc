#include "PID/DKinematicData.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DNeutralParticleHypothesis.h"

//Declare thread_local resource pools
thread_local shared_ptr<DResourcePool<DKinematicData::DKinematicInfo>> DKinematicData::dResourcePool_KinematicInfo = std::make_shared<DResourcePool<DKinematicData::DKinematicInfo>>();
thread_local shared_ptr<DResourcePool<DChargedTrackHypothesis::DTimingInfo>> DChargedTrackHypothesis::dResourcePool_TimingInfo = std::make_shared<DResourcePool<DChargedTrackHypothesis::DTimingInfo>>();
thread_local shared_ptr<DResourcePool<DChargedTrackHypothesis::DTrackingInfo>> DChargedTrackHypothesis::dResourcePool_TrackingInfo = std::make_shared<DResourcePool<DChargedTrackHypothesis::DTrackingInfo>>();
thread_local shared_ptr<DResourcePool<DNeutralParticleHypothesis::DTimingInfo>> DNeutralParticleHypothesis::dResourcePool_TimingInfo = std::make_shared<DResourcePool<DNeutralParticleHypothesis::DTimingInfo>>();

