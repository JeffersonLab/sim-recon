// $Id$
//
//    File: DNeutralParticleHypothesis.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticleHypothesis_
#define _DNeutralParticleHypothesis_

#include <vector>
#include <JANA/JObject.h>
#include <PID/DKinematicData.h>
#include <particleType.h>

using namespace std;

class DNeutralParticleHypothesis : public DKinematicData {
	public:
		JOBJECT_PUBLIC(DNeutralParticleHypothesis);

		Particle_t dPID;
		float dChiSq;
		unsigned int dNDF;
		float dFOM;

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "PID", "%d", int(dPID));
			DKinematicData::toStrings(items);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

};

#endif // _DNeutralParticleHypothesis_

