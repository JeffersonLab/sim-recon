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

		//THESE RETURN NULL IF NO MATCH TO THAT SYSTEM
		const DSCHitMatchParams* Get_SCHitMatchParams(void) const;
		const DTOFHitMatchParams* Get_TOFHitMatchParams(void) const;
		const DBCALShowerMatchParams* Get_BCALShowerMatchParams(void) const;
		const DFCALShowerMatchParams* Get_FCALShowerMatchParams(void) const;

		void Set_SCHitMatchParams(const DSCHitMatchParams& locSCHitMatchParams);
		void Set_TOFHitMatchParams(const DTOFHitMatchParams& locTOFHitMatchParams);
		void Set_BCALShowerMatchParams(const DBCALShowerMatchParams& locBCALShowerMatchParams);
		void Set_FCALShowerMatchParams(const DFCALShowerMatchParams& locFCALShowerMatchParams);

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

	private:
		//IGNORE DATA IF THE dTrack MEMBER OF THESE IS NULL: NO MATCH
		DSCHitMatchParams dSCHitMatchParams;
		DTOFHitMatchParams dTOFHitMatchParams;
		DBCALShowerMatchParams dBCALShowerMatchParams;
		DFCALShowerMatchParams dFCALShowerMatchParams;
};

//THESE RETURN NULL IF NO MATCH TO THAT SYSTEM
inline const DSCHitMatchParams* DChargedTrackHypothesis::Get_SCHitMatchParams(void) const
{
	return ((dSCHitMatchParams.dSCHit == NULL) ? NULL : &dSCHitMatchParams);
}

inline const DTOFHitMatchParams* DChargedTrackHypothesis::Get_TOFHitMatchParams(void) const
{
	return ((dTOFHitMatchParams.dTOFPoint == NULL) ? NULL : &dTOFHitMatchParams);
}

inline const DBCALShowerMatchParams* DChargedTrackHypothesis::Get_BCALShowerMatchParams(void) const
{
	return ((dBCALShowerMatchParams.dBCALShower == NULL) ? NULL : &dBCALShowerMatchParams);
}

inline const DFCALShowerMatchParams* DChargedTrackHypothesis::Get_FCALShowerMatchParams(void) const
{
	return ((dFCALShowerMatchParams.dFCALShower == NULL) ? NULL : &dFCALShowerMatchParams);
}

inline void DChargedTrackHypothesis::Set_SCHitMatchParams(const DSCHitMatchParams& locSCHitMatchParams)
{
	dSCHitMatchParams = locSCHitMatchParams;
}

inline void DChargedTrackHypothesis::Set_TOFHitMatchParams(const DTOFHitMatchParams& locTOFHitMatchParams)
{
	dTOFHitMatchParams = locTOFHitMatchParams;
}

inline void DChargedTrackHypothesis::Set_BCALShowerMatchParams(const DBCALShowerMatchParams& locBCALShowerMatchParams)
{
	dBCALShowerMatchParams = locBCALShowerMatchParams;
}

inline void DChargedTrackHypothesis::Set_FCALShowerMatchParams(const DFCALShowerMatchParams& locFCALShowerMatchParams)
{
	dFCALShowerMatchParams = locFCALShowerMatchParams;
}

#endif // _DChargedTrackHypothesis_

