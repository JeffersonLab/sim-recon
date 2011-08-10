// $Id$
//
//    File: DChargedTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_
#define _DChargedTrackHypothesis_

#include <vector>
#include <JANA/JObject.h>
#include <TRACKING/DTrackTimeBased.h>
#include <particleType.h>

using namespace std;

class DChargedTrackHypothesis : public jana::JObject {
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		const DTrackTimeBased* dTrackTimeBased;
		Particle_t dPID;
		float dProjectedTime; //Time at the track position in the DTrackTimeBased object, calculated from matching to either the FCAL, BCAL, or TOF
		float dPathLength; //Path length from the track position in the DTrackTimeBased object to the matched hit in either the FCAL, BCAL, or TOF
		float dFlightTime; //The amount of time that the track took to traverse the dPathLength
		float dChiSq;
		unsigned int dNDF;
		float dFOM;

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "PID", "%d", int(dPID));
			dTrackTimeBased->toStrings(items);
			AddString(items, "T_Proj", "%3.2f", dProjectedTime);
			AddString(items, "Path", "%3.2f", dPathLength);
			AddString(items, "TOF", "%3.2f", dFlightTime);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

};

#endif // _DChargedTrackHypothesis_

