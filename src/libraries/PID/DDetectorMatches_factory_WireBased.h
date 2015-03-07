// $Id$
//
//    File: DDetectorMatches_factory_WireBased.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DDetectorMatches_factory_WireBased_
#define _DDetectorMatches_factory_WireBased_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <PID/DDetectorMatches.h>
#include <TRACKING/DTrackWireBased.h>
#include <PID/DParticleID.h>
#include <TOF/DTOFPoint.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TMath.h>

using namespace std;
using namespace jana;

class DDetectorMatches_factory_WireBased : public jana::JFactory<DDetectorMatches>
{
	public:
		DDetectorMatches_factory_WireBased(){};
		virtual ~DDetectorMatches_factory_WireBased(){};
		const char* Tag(void){return "WireBased";}

		//called by DDetectorMatches tag=Combo factory
		DDetectorMatches* Create_DDetectorMatches(jana::JEventLoop* locEventLoop, vector<const DTrackWireBased*>& locTrackWireBasedVector);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		//matching tracks to hits/showers routines
		void MatchToTOF(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DTOFPoint*>& locTOFPoints, DDetectorMatches* locDetectorMatches) const;
		void MatchToBCAL(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DBCALShower*>& locBCALShowers, DDetectorMatches* locDetectorMatches) const;
		void MatchToFCAL(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DFCALShower*>& locFCALShowers, DDetectorMatches* locDetectorMatches) const;
		void MatchToSC(const DParticleID* locParticleID, const DTrackWireBased* locTrackWireBased, const vector<const DSCHit*>& locSCHits, DDetectorMatches* locDetectorMatches) const;

		//matching showers to tracks routines
		void MatchToTrack(const DParticleID* locParticleID, const DBCALShower* locBCALShower, const vector<const DTrackWireBased*>& locTrackWireBasedVector, DDetectorMatches* locDetectorMatches) const;
		void MatchToTrack(const DParticleID* locParticleID, const DFCALShower* locFCALShower, const vector<const DTrackWireBased*>& locTrackWireBasedVector, DDetectorMatches* locDetectorMatches) const;
};

#endif // _DDetectorMatches_factory_WireBased_
