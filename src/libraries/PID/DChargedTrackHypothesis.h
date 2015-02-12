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
#include <PID/DDetectorMatches.h>

using namespace std;

class DChargedTrackHypothesis : public DKinematicData
{
	public:
		JOBJECT_PUBLIC(DChargedTrackHypothesis);

		oid_t candidateid;   ///< id of DTrackCandidate corresponding to this track

		unsigned int dNDF_Track;
		double dChiSq_Track;
		
		unsigned int dNDF_Timing;
		double dChiSq_Timing;

		unsigned int dNDF_DCdEdx;
		double dChiSq_DCdEdx;
		
		unsigned int dNDF; //total NDF used for PID determination
		double dChiSq; //total chi-squared used for PID determination
		double dFOM; //overall FOM for PID determination

		//selected matches to detector systems
			//IGNORE DATA IF THE dTrackTimeBased MEMBER OF THESE IS NULL: NO MATCH
		DSCHitMatchParams dSCHitMatchParams;
		DTOFHitMatchParams dTOFHitMatchParams;
		DBCALShowerMatchParams dBCALShowerMatchParams;
		DFCALShowerMatchParams dFCALShowerMatchParams;

		void toStrings(vector<pair<string,string> > &items) const
		{
			AddString(items, "candidate","%d",candidateid);
			DKinematicData::toStrings(items);
			AddString(items, "Track_ChiSq", "%f", dChiSq_Track);
			AddString(items, "dEdx_ChiSq", "%f", dChiSq_DCdEdx);
			AddString(items, "TOF_ChiSq", "%f", dChiSq_Timing);
			AddString(items, "PID_ChiSq", "%f", dChiSq);
			AddString(items, "PID_FOM", "%f", dFOM);
		}
};

#endif // _DChargedTrackHypothesis_

