// $Id$
//
//    File: DTrackCandidate_factory.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_
#define _DTrackCandidate_factory_


#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;
#include <TH2F.h>
#include <TH1F.h>
#include "DTrackCandidate.h"
#include <DVector3.h>
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCSegment.h"
#include "DHelicalFit.h"
#include "DMagneticFieldStepper.h"

class DMagneticFieldMap;

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="CORE.png" width="100">
///	</A>
/// \endhtmlonly

/// Form complete list of DTrackCandidate objects using the lists formed
/// from the CDC and FDCCathodes candidate factories (DTrackCandidate_factory_CDC
/// and DTrackCandidate_factory_FDCCathodes). 
///
/// Track finding starts by looking for candidates independently in the CDC
/// and FDC. The results of those first passes are used as input here where
/// a single list is made containijng all candidates.
///
/// This will attempt to identify any candidates that should be merged into a
/// single candidate, mainly if a both a CDC and FDC candidate were found for
/// the same track.
///
/// In addition, stray CDC hits that did not belong to any candidate are
/// merged into existing candidates if possible.

class DTrackCandidate_factory:public JFactory<DTrackCandidate>{
 public:
  DTrackCandidate_factory(){
    DEBUG_HISTS=false;
    //DEBUG_HISTS=true;
  };
  ~DTrackCandidate_factory(){};
   
 protected:
  virtual jerror_t init(void);
  virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
  virtual jerror_t brun(JEventLoop* eventLoop,int runnumber);
  virtual jerror_t erun(void);
  virtual jerror_t fini(void);

  double DocaToHelix(const DCDCTrackHit *hit,double q,const DVector3 &pos,
		     const DVector3 &mom);
  double GetSenseOfRotation(DHelicalFit &fit,const DFDCPseudo *fdchit,const DVector3 &pos);
  jerror_t DoRefit(DHelicalFit &fit,
		   vector<const DFDCSegment *>segments,
		   vector<const DCDCTrackHit *>cdchits,
		   double &Bz);
  void ProjectHelixToZ(const double z,const double q,const DVector3 &mom,
		       DVector3 &pos);

  jerror_t GetPositionAndMomentum(const DFDCSegment *segment,
				  DVector3 &pos, DVector3 &mom);
  jerror_t GetPositionAndMomentum(DHelicalFit &fit,double Bz,
				  const DVector3 &origin,
				  DVector3 &pos,
				  DVector3 &mom);
  jerror_t GetPositionAndMomentum(DHelicalFit &fit,double Bz,DVector3 &pos,
				  DVector3 &mom);
  jerror_t GetPositionAndMomentum(double z,DHelicalFit &fit,
				  double Bz,DVector3 &pos,DVector3 &mom);

  void UpdatePositionAndMomentum(DTrackCandidate *can,const DFDCPseudo *fdchit,
				 DHelicalFit &fit,double Bz_avg,int axial_id);

  // Various methods for matching CDC and FDC candidates
  bool MatchMethod1(const DTrackCandidate *fdccan,
		    vector<unsigned int> &cdc_forward_ids,
		    vector<DVector3>&cdc_endplate_projections,
		    vector<unsigned int>&used_cdc_hits
		    );
  bool MatchMethod2(const DTrackCandidate *fdccan,
		    vector<unsigned int> &cdc_forward_ids,
		    vector<unsigned int>&used_cdc_hits
		    );
  bool MatchMethod3(const DTrackCandidate *cdccan,vector<int> &forward_matches,
		    vector<unsigned int>&used_cdc_hits
		    );  
  bool MatchMethod4(const DTrackCandidate *srccan,vector<int> &forward_matches,
		    int &num_fdc_cands_remaining);
  bool MatchMethod5(DTrackCandidate *can,  
		    vector<const DCDCTrackHit *>&cdchits,
		    vector<int> &forward_matches);
  void MatchMethod6(DTrackCandidate *can, 
		    vector<const DFDCSegment *>&segments,
		    vector<unsigned int>&used_cdc_hits,  
		    unsigned int &num_unmatched_cdcs
		    );
  bool MatchMethod7(DTrackCandidate *srccan,vector<int> &forward_matches,
		    int &num_fdc_cands_remaining);
  bool MatchMethod8(const DTrackCandidate *cdccan,vector<int> &forward_matches,
		    vector<unsigned int>&used_cdc_hits);
  bool MatchMethod9(unsigned int src_index,const DTrackCandidate *srccan, 
		    const DFDCSegment *segment,
		    vector<const DTrackCandidate*>&cands,
		    vector<int> &forward_matches);
  bool MatchMethod10(unsigned int src_index,const DTrackCandidate *srccan, 
		     const DFDCSegment *segment,
		     vector<const DTrackCandidate*>&cands,
		     vector<int> &forward_matches);
  bool MatchMethod11(double q,DVector3 &mypos,DVector3 &mymom,
		     DHelicalFit &fit2,const DFDCSegment *segment1,
		     const DFDCSegment *segment2);
 
 private:
  const DMagneticFieldMap *bfield;
  DMagneticFieldStepper *stepper;

  vector<const DTrackCandidate*>cdctrackcandidates;
  vector<const DTrackCandidate*>fdctrackcandidates; 
  vector<const DCDCTrackHit*>mycdchits;
  vector<DTrackCandidate *>trackcandidates;

  int DEBUG_LEVEL,MIN_NUM_HITS;
  bool DEBUG_HISTS;
  TH2F *match_dist,*match_dist_vs_p;
  TH2F *match_center_dist2;

  double FactorForSenseOfRotation;
  DVector3 cdc_endplate;
  double endplate_rmax;
  double TARGET_Z;
  int MAX_NUM_TRACK_CANDIDATES; //used to avoid memory spikes: if this # is exceeded, delete all tracks //to disable, set = -1!!
};

#endif // _DTrackCandidate_factory_

