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
#include "DTrackCandidate.h"
#include <DVector3.h>
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCSegment.h"
#include "DHelicalFit.h"

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
  
  jerror_t GetPositionAndMomentum(const DFDCSegment *segment,
				  DVector3 &pos, DVector3 &mom);
  jerror_t GetPositionAndMomentum(DHelicalFit &fit,double Bz,
				  const DVector3 &origin,
				  DVector3 &pos,
				  DVector3 &mom);
  
 protected:
  virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
  virtual jerror_t brun(JEventLoop* eventLoop,int runnumber);

  double DocaToHelix(const DCDCTrackHit *hit,double q,const DVector3 &pos,
		     const DVector3 &mom);
  double GetCharge(DHelicalFit &fit,const DFDCPseudo *fdchit,const DVector3 &pos);
  
 private:
  const DMagneticFieldMap *bfield;

  bool DEBUG_HISTS;
  TH2F *match_dist,*match_dist_vs_p;

  DVector3 cdc_endplate;
  double endplate_rmax;
  double TARGET_Z;
};

#endif // _DTrackCandidate_factory_

