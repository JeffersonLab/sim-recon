// $Id$
//
//    File: DEventRFBunch_factory_Thrown.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DEventRFBunch_factory_Thrown_
#define _DEventRFBunch_factory_Thrown_

#include <iostream>
#include <iomanip>

#include <JANA/JFactory.h>
#include <PID/DEventRFBunch.h>
#include <PID/DBeamPhoton.h>
#include <RF/DRFTime.h>
#include <TRACKING/DMCThrown.h>

using namespace std;
using namespace jana;

class DEventRFBunch_factory_Thrown : public jana::JFactory<DEventRFBunch>
{
	public:
		DEventRFBunch_factory_Thrown(){};
		~DEventRFBunch_factory_Thrown(){};
		const char* Tag(void){return "Thrown";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventRFBunch_factory_Thrown_

