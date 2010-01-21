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
#include "FDC/DFDCSegment.h"

class DMagneticFieldMap;

class DTrackCandidate_factory:public JFactory<DTrackCandidate>{
 public:
  DTrackCandidate_factory(){
    DEBUG_HISTS=false;
    //DEBUG_HISTS=true;
  };
  ~DTrackCandidate_factory(){};
  
  jerror_t GetPositionAndMomentum(const DFDCSegment *segment,
				  DVector3 &pos, DVector3 &mom);
  
 protected:
  virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
  virtual jerror_t brun(JEventLoop* eventLoop,int runnumber);
  
 private:
  const DMagneticFieldMap *bfield;

  bool DEBUG_HISTS;
  TH2F *match_dist,*match_dist_vs_p;

  DVector3 cdc_endplate;
  double endplate_rmax;
};

#endif // _DTrackCandidate_factory_

