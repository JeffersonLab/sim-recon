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
#include <TRACKING/DTrackingData.h>
#include <TRACKING/DTrackFitter.h>

class DReferenceTrajectory;


class DTrackWireBased:public DTrackingData{
	public:
		JOBJECT_PUBLIC(DTrackWireBased);
		
		oid_t candidateid;	///< which DTrackCandidate this came from
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		vector<DTrackFitter::pull_t> pulls;	///< Holds pulls used in chisq calc. (not including off-diagonals)
		double FOM; //confidence level

		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track

      bool IsSmoothed; // Boolean value to indicate whether the smoother was run succesfully over this track.

		// Hit CDC Rings & FDC Planes
		// use the DParticleID Get_CDCRings & Get_FDCPlanes functions to extract the information from these
		unsigned int dCDCRings; //CDC rings where the track has an associated DCDCTrackHit //rings correspond to bits (1 -> 28)
		unsigned int dFDCPlanes; //FDC planes where the track has an associated DFDCPseudoHit //planes correspond to bits (1 -> 24)

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "candidate", "%d", candidateid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DTrackWireBased_

