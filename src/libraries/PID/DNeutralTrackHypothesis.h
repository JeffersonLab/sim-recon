// $Id$
//
//    File: DNeutralTrackHypothesis.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralTrackHypothesis_
#define _DNeutralTrackHypothesis_

#include <vector>
#include <JANA/JObject.h>
#include <PID/DKinematicData.h>
#include <particleType.h>

using namespace std;

class DNeutralTrackHypothesis : public jana::JObject {
	public:
		JOBJECT_PUBLIC(DNeutralTrackHypothesis);

		DKinematicData* dKinematicData;
		Particle_t dPID;
		float dProjectedTime; //Time at the track position in the DKinematicData object, which is the time of the associated DVertex
		float dPathLength; //Path length from the vertex to the matched hit in either the FCAL or BCAL
		float dFlightTime; //The amount of time that the track takes to traverse the dPathLength
		float dChiSq;
		unsigned int dNDF;
		float dFOM;
};

#endif // _DNeutralTrackHypothesis_

