// $Id$
//
//    File: DChargedTrack.cc
// Created: Tue Mar  6 14:29:24 EST 2012
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "PID/DChargedTrack.h"

int DChargedTrack::Get_Charge(void) const
{
	const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_BestFOM();
	return ((locChargedTrackHypothesis == NULL) ? 0 : ParticleCharge(locChargedTrackHypothesis->PID()));
}

const DChargedTrackHypothesis* DChargedTrack::Get_BestFOM(void) const
{
	if(dChargedTrackHypotheses.empty())
		return NULL;
	double locBestFOM = -2.0;
	const DChargedTrackHypothesis* locBestChargedTrackHypothesis = dChargedTrackHypotheses[0];
	for(size_t loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i)
	{
		if(dChargedTrackHypotheses[loc_i]->dFOM > locBestFOM)
		{
			locBestChargedTrackHypothesis = dChargedTrackHypotheses[loc_i];
			locBestFOM = locBestChargedTrackHypothesis->dFOM;
		}
	}
	return locBestChargedTrackHypothesis;
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

