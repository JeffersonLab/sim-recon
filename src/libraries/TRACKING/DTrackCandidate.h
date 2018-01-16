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

#include <TRACKING/DTrackingData.h>
#include <TRACKING/DTrackFitter.h>

class DReferenceTrajectory;

#define MAX_IHITS 256

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="CORE.png" width="100">
///	</A>
/// \endhtmlonly

/// DTrackCandidate objects are the result of track finding and the
/// input to track fitting. Several algorithms exist for making
/// these and then merging them. For the default, see
/// DTrackCandidate_factory .

class DTrackCandidate:public DTrackingData{
	public:
		JOBJECT_PUBLIC(DTrackCandidate);
		
		DTrackCandidate():chisq(0),Ndof(0),rt(0),IsSmoothed(0){}

		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit

		const DReferenceTrajectory *rt; ///< pointer to reference trjectory representing this track (if any)
		
		vector<DTrackFitter::pull_t> pulls; // vector of residuals and other track-related quantities 
		
		vector<int>used_cdc_indexes;
		vector<int>used_fdc_indexes;

		// Circle fit data
		double xc,yc,rc;

		// Hit CDC Rings & FDC Planes
		// use the DParticleID Get_CDCRings & Get_FDCPlanes functions to extract the information from these
		unsigned int dCDCRings; //CDC rings where the track has an associated DCDCTrackHit //rings correspond to bits (1 -> 28)
		unsigned int dFDCPlanes; //FDC planes where the track has an associated DFDCPseudoHit //planes correspond to bits (1 -> 24)

      bool IsSmoothed; // Boolean value to indicate whether the smoother was run succesfully over this track.

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "id", "0x%x", id);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
		}
};

#endif // _DTrackCandidate_

