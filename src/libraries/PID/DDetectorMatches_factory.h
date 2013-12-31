// $Id$
//
//    File: DDetectorMatches_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DDetectorMatches_factory_
#define _DDetectorMatches_factory_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <PID/DDetectorMatches.h>
#include <TRACKING/DTrackTimeBased.h>
#include <PID/DParticleID.h>
#include <TOF/DTOFPoint.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TMath.h>

using namespace std;
using namespace jana;

class DDetectorMatches_factory : public jana::JFactory<DDetectorMatches>
{
	public:
		DDetectorMatches_factory(){};
		~DDetectorMatches_factory(){};

		//called by DDetectorMatches tag=Combo factory
		DDetectorMatches* Create_DDetectorMatches(jana::JEventLoop* locEventLoop, vector<const DTrackTimeBased*>& locTrackTimeBasedVector) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		//matching tracks to hits/showers routines
		void MatchToTOF(const DTrackTimeBased* locTrackTimeBased, const vector<const DTOFPoint*>& locTOFPoints, DDetectorMatches* locDetectorMatches) const;
		void MatchToBCAL(const DTrackTimeBased* locTrackTimeBased, const vector<const DBCALShower*>& locBCALShowers, DDetectorMatches* locDetectorMatches) const;
		void MatchToFCAL(const DTrackTimeBased* locTrackTimeBased, const vector<const DFCALShower*>& locFCALShowers, DDetectorMatches* locDetectorMatches) const;
		void MatchToSC(const DTrackTimeBased* locTrackTimeBased, const vector<const DSCHit*>& locSCHits, DDetectorMatches* locDetectorMatches) const;

		//matching showers to tracks routines
		void MatchToTrack(const DBCALShower* locBCALShower, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, DDetectorMatches* locDetectorMatches) const;
		void MatchToTrack(const DFCALShower* locFCALShower, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, DDetectorMatches* locDetectorMatches) const;

		//set flight-time/p correlations
		void Set_SCFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
		void Set_TOFFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
		void Set_BCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;
		void Set_FCALFlightTimePCorrelation(const DTrackTimeBased* locTrackTimeBased, DDetectorMatches* locDetectorMatches) const;

		const DParticleID* dParticleID;
};

#endif // _DDetectorMatches_factory_

