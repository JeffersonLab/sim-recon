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

#include <TMath.h>

#include <JANA/JFactory.h>

#include <DVector3.h>
#include <DMatrixDSym.h>
#include <DMatrix.h>

#include <PID/DDetectorMatches.h>
#include <PID/DParticleID.h>
#include <PID/DEventRFBunch.h>
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

		bool Find_TimeFOMPairs_Hits(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, double> >& locTimeFOMPairs);
		bool Find_TimeFOMPairs_T0(const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, double> >& locTimeFOMPairs);
		int Find_BestRFBunchShift(double locRFHitTime, const vector<pair<double, double> >& locTimeFOMPairs);

	private:
		const DParticleID* dParticleID;

		double dRFBunchFrequency;
		DVector3 dTargetCenter;

		double dMinTrackingFOM;
		double dMinVertexZ;
		double dMaxVertexZ;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventRFBunch_factory_

