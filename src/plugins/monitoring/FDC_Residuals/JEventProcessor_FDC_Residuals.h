// $Id$
//
//    File: JEventProcessor_FDC_Residuals.cc
// Created: Wed Nov 30 10:08:00 EDT 2016
// Creator: aaustreg
//

#ifndef _JEventProcessor_FDC_Residuals_
#define _JEventProcessor_FDC_Residuals_

//#include <pthread.h>
#include <map>
#include <vector>
#include <deque>
using namespace std;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>

#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackCandidate_factory_StraightLine.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackWireBased.h>
#include <PID/DChargedTrack.h>
#include <PID/DDetectorMatches.h>
#include <FDC/DFDCHit.h>

class JEventProcessor_FDC_Residuals:public jana::JEventProcessor{
 public:
  JEventProcessor_FDC_Residuals();
  ~JEventProcessor_FDC_Residuals();
  const char* className(void){return "JEventProcessor_FDC_Residuals";}
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  DGeometry * dgeom;
  bool dIsNoFieldFlag;
  vector< vector< DFDCWire * > > fdcwires; // FDC Wires Referenced by [layer 1-24][wire 1-96]
  vector<double> fdcz; // FDC z positions
  vector<double> fdcrmin; // FDC inner radii
  double fdcrmax; // FDC outer radius
};

#endif // _JEventProcessor_FDC_Residuals_

