// $Id: DTrackWireBased.h 4927 2009-03-11 13:20:33Z staylor $
//
//    File: DTrackWireBased.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackWireBased_
#define _DTrackWireBased_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DKinematicData.h>
#include <TRACKING/DTrackFitter.h>

class DReferenceTrajectory;

typedef struct{
  double t0,t0_sigma;
  DetectorSystem_t system;
}DStartTime_t;

class DTrackWireBased:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DTrackWireBased);
		
		oid_t candidateid;	///< which DTrackCandidate this came from
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		vector<DTrackFitter::pull_t> pulls;	///< Holds pulls used in chisq calc. (not including off-diagonals)

		vector<DStartTime_t>start_times;

		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "candidate", "%d", candidateid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DTrackWireBased_

