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

class DTrackWireBased:public DTrackingData{
	public:
		JOBJECT_PUBLIC(DTrackWireBased);
		
		oid_t candidateid;	///< which DTrackCandidate this came from
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		vector<DTrackFitter::pull_t> pulls;	///< Holds pulls used in chisq calc. (not including off-diagonals)
		map<DetectorSystem_t,vector<DTrackFitter::Extrapolation_t> >extrapolations;
		
		double FOM; //confidence level

		bool GetProjection(DetectorSystem_t detector,DVector3 &pos,
				   DVector3 *mom=nullptr,double *t=nullptr) const;

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

inline bool DTrackWireBased::GetProjection(DetectorSystem_t detector,
					   DVector3 &pos,
					   DVector3 *mom,double *t) const{
 
  if (detector>SYS_BCAL && extrapolations.at(detector).size()>0){
    DTrackFitter::Extrapolation_t extrapolation=extrapolations.at(detector)[0];
    pos=extrapolation.position;
    if (mom){
      *mom=extrapolation.momentum;
    }
    if (t){
      *t=extrapolation.t;
    }
    return true;
  }
  // Set defaults that are clearly unreasonable, since the projection did not work! 
  pos.SetXYZ(0.,0.,-100.);
  if (mom){
    mom->SetXYZ(0.,0.,0.);
  }
  if (t){
    *t=-1000.;
  }  
  return false;
}

#endif // _DTrackWireBased_

