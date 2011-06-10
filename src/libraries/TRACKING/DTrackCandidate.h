// $Id$
//
//    File: DTrackCandidate.h
// Created: Sun Apr  3 12:38:16 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_
#define _DTrackCandidate_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

#include "PID/DKinematicData.h"

class DReferenceTrajectory;

#define MAX_IHITS 256

/// DTrackCandidate objects are the result of track finding and the
/// input to track fitting. Several algorithms exist for making
/// these and then merging them. For the default, see
/// DTrackCandidate_factory .

class DTrackCandidate:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DTrackCandidate);
		
		DTrackCandidate():chisq(0),Ndof(0){}

		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit

		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track (if any)

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "id", "0x%x", id);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DTrackCandidate_

