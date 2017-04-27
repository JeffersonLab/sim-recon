// $Id$
//
//    File: JEventProcessor_FDCProjectionResiduals.h
// Created: Wed Oct 26 14:07:16 EDT 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FDCProjectionResiduals_
#define _JEventProcessor_FDCProjectionResiduals_

#include <set>
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
#include <PID/DParticleID.h>
#include <PID/DDetectorMatches.h>
#include <CDC/DCDCTrackHit.h>

class JEventProcessor_FDCProjectionResiduals:public jana::JEventProcessor{
	public:
		JEventProcessor_FDCProjectionResiduals();
		~JEventProcessor_FDCProjectionResiduals();
		const char* className(void){return "JEventProcessor_FDCProjectionResiduals";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

      DGeometry * dgeom;
      bool dIsNoFieldFlag;
      double dTargetCenterZ;
      double dTargetLength;

      double MAX_DRIFT_TIME;
      int PLANE_TO_SKIP;
      double dMinTrackingFOM;

      vector< vector< DCDCWire * > > cdcwires; // CDC Wires Referenced by [ring][straw]
      vector<vector<double> >max_sag;
      vector<vector<double> >sag_phi_offset;

      vector<double> cdc_drift_table;
      double long_drift_func[3][3];
      double short_drift_func[3][3];
      double endplate_z,endplate_dz,endplate_rmin,endplate_rmax;

      bool Expect_Hit(const DTrackTimeBased* thisTimeBasedTrack, DCDCWire* wire, double distanceToWire, double& delta, double& dz);
      double CDCDriftDistance(double delta, double t);
      unsigned int Locate(vector<double>&xx,double x);
      double GetDOCA(DVector3, DVector3, DVector3, DVector3, DVector3&, DVector3&);
      
};

#endif // _JEventProcessor_FDCProjectionResiduals_

