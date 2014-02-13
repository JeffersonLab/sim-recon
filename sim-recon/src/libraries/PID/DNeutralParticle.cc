// $Id$
//
//    File: DNeutralParticle.cc
// Created: Tue Mar  6 14:29:24 EST 2012
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "PID/DNeutralParticle.h"

const DNeutralParticleHypothesis* DNeutralParticle::Get_BestFOM(void) const
{
	if(dNeutralParticleHypotheses.empty())
		return NULL;
	double locBestFOM = -2.0;
	const DNeutralParticleHypothesis* locBestNeutralParticleHypotheses = dNeutralParticleHypotheses[0];
	for(size_t loc_i = 0; loc_i < dNeutralParticleHypotheses.size(); ++loc_i)
	{
		if(dNeutralParticleHypotheses[loc_i]->dFOM > locBestFOM)
		{
			locBestNeutralParticleHypotheses = dNeutralParticleHypotheses[loc_i];
			locBestFOM = locBestNeutralParticleHypotheses->dFOM;
		}
	}
	return locBestNeutralParticleHypotheses;
}

const DNeutralParticleHypothesis* DNeutralParticle::Get_Hypothesis(Particle_t locPID) const
{
	for(unsigned int loc_i = 0; loc_i < dNeutralParticleHypotheses.size(); ++loc_i)
	{
		if(dNeutralParticleHypotheses[loc_i]->PID() == locPID)
			return dNeutralParticleHypotheses[loc_i];
	}
	return NULL;
}

