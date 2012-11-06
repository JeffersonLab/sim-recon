// $Id$
//
//    File: DChargedTrack.cc
// Created: Tue Mar  6 14:29:24 EST 2012
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "PID/DChargedTrack.h"

double DChargedTrack::Get_Charge(void) const
{
	if(dChargedTrackHypotheses.size() == 0)
		return 0.0;
	return dChargedTrackHypotheses[0]->charge();
}

const DChargedTrackHypothesis* DChargedTrack::Get_BestFOM(void) const
{
	if(dChargedTrackHypotheses.size() == 0)
		return NULL;
	return dChargedTrackHypotheses[0];
}

const DChargedTrackHypothesis* DChargedTrack::Get_Hypothesis(Particle_t locPID) const
{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		if(dChargedTrackHypotheses[loc_i]->PID() == locPID)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

