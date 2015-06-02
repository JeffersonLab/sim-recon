// $Id$
//
//    File: DEventRFBunch_factory_Calibrations.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DEventRFBunch_factory_Calibrations_
#define _DEventRFBunch_factory_Calibrations_

#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <deque>
#include <vector>

#include <TMath.h>

#include <JANA/JFactory.h>

#include <DVector3.h>
#include <DMatrixDSym.h>
#include <DMatrix.h>

#include <TTAB/DTTabUtilities.h>
#include <PID/DDetectorMatches.h>
#include <PID/DParticleID.h>
#include <PID/DEventRFBunch.h>
#include <RF/DRFTDCDigiTime.h>
#include <TRACKING/DTrackWireBased.h>
#include <START_COUNTER/DSCHit.h>

#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

#include <ANALYSIS/DCutActions.h>

using namespace std;
using namespace jana;

class DEventRFBunch_factory_Calibrations : public jana::JFactory<DEventRFBunch>
{
	public:
		DEventRFBunch_factory_Calibrations(){};
		~DEventRFBunch_factory_Calibrations(){};
		const char* Tag(void){return "Calibrations";}

	private:

		void Select_GoodTracks(JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locSelectedWireBasedTracks) const;
		jerror_t Select_RFBunch(JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locTrackWireBasedVector, double locRFTime);

		bool Find_TrackTimes_SC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackWireBased*>& locTrackWireBasedVector, vector<pair<double, const JObject*> >& locTimes) const;
		int Conduct_Vote(JEventLoop* locEventLoop, double locRFTime, vector<pair<double, const JObject*> >& locTimes, int& locHighestNumVotes);

		int Find_BestRFBunchShifts(double locRFHitTime, const vector<pair<double, const JObject*> >& locTimes, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts);
		int Break_TieVote_Tracks(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts);

		jerror_t Create_NaNRFBunch(void);

		const DParticleID* dParticleID;

		double dRFBunchPeriod;
		DVector3 dTargetCenter;

		DetectorSystem_t dRFTDCSourceSystem;
		double dMinTrackingFOM;
		unsigned int dMinHitRingsPerCDCSuperlayer;
		unsigned int dMinHitPlanesPerFDCPackage;

		DCutAction_TrackHitPattern *dCutAction_TrackHitPattern;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventRFBunch_factory_Calibrations_

