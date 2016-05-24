// $Id$
//
//    File: DNeutralParticle.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticle_
#define _DNeutralParticle_

#include <vector>
#include <JANA/JObject.h>
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DNeutralShower.h>
#include <particleType.h>

using namespace std;

class DNeutralParticle : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DNeutralParticle);

		vector<const DNeutralParticleHypothesis*> dNeutralParticleHypotheses;
		const DNeutralShower* dNeutralShower;

		const DNeutralParticleHypothesis* Get_BestFOM(void) const;
		const DNeutralParticleHypothesis* Get_Hypothesis(Particle_t locPID) const;

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "Nhypotheses", "%d", dNeutralParticleHypotheses.size());
		}
};

inline const DNeutralParticleHypothesis* DNeutralParticle::Get_BestFOM(void) const
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

inline const DNeutralParticleHypothesis* DNeutralParticle::Get_Hypothesis(Particle_t locPID) const
{
	for(unsigned int loc_i = 0; loc_i < dNeutralParticleHypotheses.size(); ++loc_i)
	{
		if(dNeutralParticleHypotheses[loc_i]->PID() == locPID)
			return dNeutralParticleHypotheses[loc_i];
	}
	return NULL;
}

#endif // _DNeutralParticle_

