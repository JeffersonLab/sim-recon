// $Id$
//
//    File: DNeutralParticle.cc
// Created: Tue Mar  6 14:29:24 EST 2012
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "PID/DNeutralParticle.h"

const DNeutralParticleHypothesis* DNeutralParticle::Get_BestFOM(void) const{
	if(dNeutralParticleHypotheses.size() == 0)
		return NULL;
	return dNeutralParticleHypotheses[0];
}

const DNeutralParticleHypothesis* DNeutralParticle::Get_BestPhoton(void) const{
	double locBestFOM = -1.0;
	const DNeutralParticleHypothesis* locBestPhotonHypothesis = NULL;
	for(unsigned int loc_i = 0; loc_i < dNeutralParticleHypotheses.size(); ++loc_i){
		if(dNeutralParticleHypotheses[loc_i]->dPID == Gamma){
			if(dNeutralParticleHypotheses[loc_i]->dFOM > locBestFOM){
				locBestPhotonHypothesis = dNeutralParticleHypotheses[loc_i];
				locBestFOM = dNeutralParticleHypotheses[loc_i]->dFOM;
			}
		}
	}
	return locBestPhotonHypothesis;
}

const DNeutralParticleHypothesis* DNeutralParticle::Get_BestNeutron(void) const{
	double locBestFOM = -1.0;
	const DNeutralParticleHypothesis* locBestNeutronHypothesis = NULL;
	for(unsigned int loc_i = 0; loc_i < dNeutralParticleHypotheses.size(); ++loc_i){
		if(dNeutralParticleHypotheses[loc_i]->dPID == Neutron){
			if(dNeutralParticleHypotheses[loc_i]->dFOM > locBestFOM){
				locBestNeutronHypothesis = dNeutralParticleHypotheses[loc_i];
				locBestFOM = dNeutralParticleHypotheses[loc_i]->dFOM;
			}
		}
	}
	return locBestNeutronHypothesis;
}


