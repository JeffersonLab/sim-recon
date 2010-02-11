// $Id$
//
//    File: DTrackTimeBased.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_
#define _DTrackTimeBased_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <PID/DKinematicData.h>

class DReferenceTrajectory;


class DTrackTimeBased:public DKinematicData{
	public:
		JOBJECT_PUBLIC(DTrackTimeBased);
		
		oid_t trackid;			///< id of DTrack this came from
		oid_t candidateid;   /// < id of DTrackCandidate corresponding to this track
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		//float dE;				///< Total energy deposited in straws
		//float ds;				///< Total pathlength through straws contributing to dE
		//float err_dE;			///< Error on value of dE
		//float err_ds;			///< Error on value of ds
		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track

		double FOM;
		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "trackid", "0x%x", trackid);
			AddString(items, "candidateid","%d",candidateid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
			AddString(items, "FOM", "%f",FOM);
		}
};

#endif // _DTrackTimeBased_

