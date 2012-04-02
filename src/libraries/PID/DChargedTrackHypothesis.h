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

		oid_t candidateid;   ///< id of DTrackCandidate corresponding to this track

		Particle_t dPID;

		unsigned int dNDF_Track;
		float dChiSq_Track;
		
		unsigned int dNDF_Timing;
		float dChiSq_Timing;

		unsigned int dNDF_DCdEdx;
		float dChiSq_DCdEdx;

		unsigned int dNDF; //total NDF used for PID determination
		float dChiSq; //total chi-squared used for PID determination
		float dFOM; //overall FOM for PID determination

		void toStrings(vector<pair<string,string> > &items) const{
			AddString(items, "PID", "%d", int(dPID));
			DKinematicData::toStrings(items);	
			AddString(items, "Track_ChiSq", "%f", dChiSq_Track);
			AddString(items, "dEdx_ChiSq", "%f", dChiSq_DCdEdx);
			AddString(items, "TOF_ChiSq", "%f", dChiSq_Timing);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}

};

#endif // _DChargedTrackHypothesis_

