// $Id$
//
//    File: JEventProcessor_FDC_Efficiency.cc
// Created: Thu May 26 15:57:38 EDT 2016
// Creator: aaustreg
//

#ifndef _JEventProcessor_FDC_Efficiency_
#define _JEventProcessor_FDC_Efficiency_

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

class JEventProcessor_FDC_Efficiency:public jana::JEventProcessor{
 public:
  JEventProcessor_FDC_Efficiency();
  ~JEventProcessor_FDC_Efficiency();
  const char* className(void){return "JEventProcessor_FDC_Efficiency";}
  
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

#endif // _JEventProcessor_FDC_Efficiency_

