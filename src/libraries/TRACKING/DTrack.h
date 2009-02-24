// $Id$
//
//    File: DTrack.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_
#define _DTrack_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"
#include "PID/DKinematicData.h"

class DReferenceTrajectory;

class DTrack:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DTrack);
		
		oid_t candidateid;	///< id of DTrackCandidate this came from
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit

		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "candidateid", "0x%x", candidateid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DTrack_

