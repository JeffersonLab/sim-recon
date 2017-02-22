// $Id$
//
//    File: DChargedTrack_factory_PreSelect.h
// Created: Mon Dec  7 14:29:24 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrack_factory_PreSelect_
#define _DChargedTrack_factory_PreSelect_

#include <JANA/JFactory.h>
#include <PID/DChargedTrack.h>
#include <PID/DChargedTrackHypothesis.h>
#include <TRACKING/DTrackTimeBased.h>

using namespace std;
using namespace jana;

class DChargedTrack_factory_PreSelect : public jana::JFactory<DChargedTrack>
{
	public:
		DChargedTrack_factory_PreSelect(){};
		~DChargedTrack_factory_PreSelect(){};
		const char* Tag(void){return "PreSelect";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		bool Cut_HasDetectorMatch(const DChargedTrackHypothesis* locChargedTrackHypothesis, const DDetectorMatches* locDetectorMatches) const;
		bool Cut_TrackingFOM(const DChargedTrackHypothesis* locChargedTrackHypothesis) const;

		//Command-line values will override these
		double dMinTrackingFOM; //PRESELECT:MIN_TRACKING_FOM 
		bool dHasDetectorMatchFlag; //PRESELECT:HAS_DETECTOR_MATCH_FLAG: if true/false, do/don't require tracks to have a detector match
};

#endif // _DChargedTrack_factory_PreSelect_

