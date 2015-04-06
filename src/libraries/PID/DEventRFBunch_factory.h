// $Id$
//
//    File: DEventRFBunch_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DEventRFBunch_factory_
#define _DEventRFBunch_factory_

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

#include <PID/DDetectorMatches.h>
#include <PID/DParticleID.h>
#include <PID/DEventRFBunch.h>
#include <PID/DNeutralShower.h>
#include <PID/DRFTime.h>
#include <TRACKING/DTrackTimeBased.h>
#include <START_COUNTER/DSCHit.h>
#include <BCAL/DBCALShower.h>
#include <TOF/DTOFPoint.h>
#include <PID/DBeamPhoton.h>

#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

using namespace std;
using namespace jana;

class DEventRFBunch_factory : public jana::JFactory<DEventRFBunch>
{
	public:
		DEventRFBunch_factory(){};
		~DEventRFBunch_factory(){};

		bool Get_RFTimeGuess(JEventLoop* locEventLoop, double& locRFTimeGuess, double& locRFVariance, DetectorSystem_t& locTimeSource) const;

	private:

		void Select_GoodTracks(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locSelectedTimeBasedTracks) const;
		jerror_t Select_RFBunch(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector, const DRFTime* locRFTime);
		int Conduct_Vote(JEventLoop* locEventLoop, double locRFTime, vector<pair<double, const JObject*> >& locTimes, bool locUsedTracksFlag, int& locHighestNumVotes);

		bool Find_TrackTimes_SC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, const JObject*> >& locTimes) const;
		bool Find_TrackTimes_NonSC(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, const JObject*> >& locTimes);
		bool Find_NeutralTimes(JEventLoop* locEventLoop, vector<pair<double, const JObject*> >& locTimes);

		int Find_BestRFBunchShifts(double locRFHitTime, const vector<pair<double, const JObject*> >& locTimes, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts);

		bool Break_TieVote_BeamPhotons(vector<const DBeamPhoton*>& locBeamPhotons, double locRFTime, map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts, int locHighestNumVotes);
		int Break_TieVote_Tracks(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts);
		int Break_TieVote_Neutrals(map<int, vector<const JObject*> >& locNumRFBucketsShiftedMap, set<int>& locBestRFBunchShifts);

		jerror_t Select_RFBunch_NoRFTime(JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector);

		void Get_RFTimeGuess(vector<pair<double, const JObject*> >& locTimes, double& locRFTimeGuess, double& locRFVariance) const;

		jerror_t Create_NaNRFBunch(void);

		const DParticleID* dParticleID;

		double dRFBunchPeriod;
		DVector3 dTargetCenter;

		double dMinTrackingFOM;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventRFBunch_factory_

