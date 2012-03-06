// $Id$
//
//    File: DChargedTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_
#define _DChargedTrackHypothesis_

#include <vector>
#include <PID/DKinematicData.h>
#include <particleType.h>

using namespace std;

class DChargedTrackHypothesis : public DKinematicData {
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		Particle_t dPID;

		unsigned int dNDF_Timing;
		float dChiSq_Timing;

		unsigned int dNDF_DCdEdx;
		float dChiSq_DCdEdx;

		unsigned int dNDF;
		float dChiSq;
		float dFOM;

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "PID", "%d", int(dPID));
			DKinematicData::toStrings(items);
			AddString(items, "dEdx_ChiSq", "%f", dChiSq_DCdEdx);
			AddString(items, "TOF_ChiSq", "%f", dChiSq_Timing);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

};

#endif // _DChargedTrackHypothesis_

