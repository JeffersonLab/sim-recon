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
#include "DTrackingData.h"
#include "DTrackFitter.h"
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"

using namespace jana;
using namespace std;

class DTrackTimeBased:public DTrackingData{
	public:
		JOBJECT_PUBLIC(DTrackTimeBased);
		
		double dEdx(void) const{return ((dNumHitsUsedFordEdx_CDC >= dNumHitsUsedFordEdx_FDC) ? ddEdx_CDC_amp : ddEdx_FDC);}
		typedef struct{
		  unsigned int inner_layer;
		  unsigned int outer_layer;
		  unsigned int total_hits;
		}hit_usage_t;

		hit_usage_t cdc_hit_usage;
		hit_usage_t fdc_hit_usage;

		unsigned int measured_fdc_hits_on_track;
		unsigned int measured_cdc_hits_on_track;
		unsigned int potential_fdc_hits_on_track;
		unsigned int potential_cdc_hits_on_track;

		oid_t trackid;			///< id of DTrackWireBased corresponding to this track
		oid_t candidateid;   ///< id of DTrackCandidate corresponding to this track
		float chisq;			///< Chi-squared for the track (not chisq/dof!)
		int Ndof;				///< Number of degrees of freedom in the fit
		vector<DTrackFitter::pull_t> pulls;	///< Holds pulls used in chisq calc. (not including off-diagonals)
		map<DetectorSystem_t,vector<DTrackFitter::Extrapolation_t> >extrapolations;
		int flags;
		enum DTrackTimeBased_flag_t{
		  FLAG__GOODFIT=0,
		  FLAG__USED_WIREBASED_FIT=1,
		  FLAG__USED_OTHER_HYPOTHESIS=2,
		};
		

		bool GetProjection(DetectorSystem_t detector,DVector3 &pos,
				   DVector3 *mom=nullptr,double *t=nullptr) const;


      bool IsSmoothed; // Boolean value to indicate whether the smoother was run succesfully over this track.

		typedef struct{
		  double t0,t0_sigma;
		  DetectorSystem_t system;
		}DStartTime_t;
		vector<DStartTime_t>start_times;

		double FOM;

		double ddEdx_FDC;
		double ddx_FDC;
		unsigned int dNumHitsUsedFordEdx_FDC;
		double ddEdx_CDC,ddEdx_CDC_amp;
		double ddx_CDC,ddx_CDC_amp;
		unsigned int dNumHitsUsedFordEdx_CDC;

		// Hit CDC Rings & FDC Planes
		// use the DParticleID Get_CDCRings & Get_FDCPlanes functions to extract the information from these
		unsigned int dCDCRings; //CDC rings where the track has an associated DCDCTrackHit //rings correspond to bits (1 -> 28)
		unsigned int dFDCPlanes; //FDC planes where the track has an associated DFDCPseudoHit //planes correspond to bits (1 -> 24)

		// Matching to MC: Highest % of track hits matched to a thrown
		int dMCThrownMatchMyID; //MC track match myid (-1 if somehow no match)
		int dNumHitsMatchedToThrown;

		void toStrings(vector<pair<string,string> > &items)const{
			DKinematicData::toStrings(items);
			AddString(items, "candidate","%d",candidateid);
			//AddString(items, "wirebased","%d",trackid);
			AddString(items, "chisq", "%f", chisq);
			AddString(items, "Ndof", "%d", Ndof);
			AddString(items, "FOM", "%f",(float)FOM);
			AddString(items, "Flags","%d",flags);
			//AddString(items, "MCMatchID", "%d",dMCThrownMatchMyID);
			//AddString(items, "#HitsMCMatched", "%d",dNumHitsMatchedToThrown);
		}
};

size_t Get_NumTrackHits(const DTrackTimeBased* locTrackTimeBased);
inline size_t Get_NumTrackHits(const DTrackTimeBased* locTrackTimeBased)
{
	vector<const DCDCTrackHit*> locCDCHits;
	locTrackTimeBased->Get(locCDCHits);
	vector<const DFDCPseudo*> locFDCHits;
	locTrackTimeBased->Get(locFDCHits);

	size_t locNumHits = locCDCHits.size() + locFDCHits.size();
	if(locNumHits > 0)
		return locNumHits;

	return locTrackTimeBased->Ndof + 5; //is WRONG because FDC DoF != FDC Hits
}

inline bool DTrackTimeBased::GetProjection(DetectorSystem_t detector,
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

  
  return false;
}

#endif // _DTrackTimeBased_

