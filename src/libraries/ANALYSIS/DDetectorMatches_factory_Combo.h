// $Id$
//
//    File: DDetectorMatches_factory_Combo.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DDetectorMatches_factory_Combo_
#define _DDetectorMatches_factory_Combo_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <PID/DDetectorMatches_factory.h>

#include <PID/DDetectorMatches.h>
#include <TRACKING/DTrackTimeBased.h>
#include <PID/DParticleID.h>
#include <TOF/DTOFPoint.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>

using namespace std;
using namespace jana;

class DDetectorMatches_factory_Combo : public jana::JFactory<DDetectorMatches>
{
	public:
		DDetectorMatches_factory_Combo(){use_factory = 1;}; //prevents JANA from searching the input file for these objects
		~DDetectorMatches_factory_Combo(){};
		const char* Tag(void){return "Combo";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DDetectorMatches_factory* dDetectorMatchesFactory;
};

#endif // _DDetectorMatches_factory_Combo_

