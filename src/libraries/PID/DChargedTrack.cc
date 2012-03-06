// $Id$
//
//    File: DChargedTrack.cc
// Created: Tue Mar  6 14:29:24 EST 2012
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "PID/DChargedTrack.h"

double DChargedTrack::Get_Charge(void) const{
	if(dChargedTrackHypotheses.size() == 0)
		return 0.0;
	return dChargedTrackHypotheses[0]->charge();
}

const DChargedTrackHypothesis* DChargedTrack::Get_BestFOM(void) const{
	if(dChargedTrackHypotheses.size() == 0)
		return NULL;
	return dChargedTrackHypotheses[0];
}

const DChargedTrackHypothesis* DChargedTrack::Get_Proton(void) const{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i){
		if(dChargedTrackHypotheses[loc_i]->dPID == Proton)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

const DChargedTrackHypothesis* DChargedTrack::Get_PiPlus(void) const{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i){
		if(dChargedTrackHypotheses[loc_i]->dPID == PiPlus)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

const DChargedTrackHypothesis* DChargedTrack::Get_KPlus(void) const{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i){
		if(dChargedTrackHypotheses[loc_i]->dPID == KPlus)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

const DChargedTrackHypothesis* DChargedTrack::Get_PiMinus(void) const{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i){
		if(dChargedTrackHypotheses[loc_i]->dPID == PiMinus)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

const DChargedTrackHypothesis* DChargedTrack::Get_KMinus(void) const{
	for(unsigned int loc_i = 0; loc_i < dChargedTrackHypotheses.size(); ++loc_i){
		if(dChargedTrackHypotheses[loc_i]->dPID == KMinus)
			return dChargedTrackHypotheses[loc_i];
	}
	return NULL;
}

