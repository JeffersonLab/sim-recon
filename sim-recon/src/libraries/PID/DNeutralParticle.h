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
#include <particleType.h>

using namespace std;

class DNeutralParticle:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DNeutralParticle);

		vector<const DNeutralParticleHypothesis *> dNeutralParticleHypotheses;

		const DNeutralParticleHypothesis* Get_BestFOM(void) const;
		const DNeutralParticleHypothesis* Get_BestPhoton(void) const;
		const DNeutralParticleHypothesis* Get_BestNeutron(void) const;

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "Nhypotheses", "%d", dNeutralParticleHypotheses.size());
		}

};

#endif // _DNeutralParticle_

