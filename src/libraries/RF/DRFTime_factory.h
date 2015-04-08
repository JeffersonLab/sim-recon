// $Id$
//
//    File: DRFTime_factory.h
// Created: Mon Mar 30 10:51:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DRFTime_factory_
#define _DRFTime_factory_

#include <limits>
#include <iostream>

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DRFTime.h"
#include "DRFDigiTime.h"
#include "DRFTDCDigiTime.h"

using namespace std;
using namespace jana;

class DRFTime_factory : public jana::JFactory<DRFTime>
{
	public:
		DRFTime_factory(){};
		~DRFTime_factory(){};

		double Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo) const;

		double Convert_TDCToTime(const DRFTDCDigiTime* locRFTDCDigiTime, const DTTabUtilities* locTTabUtilities) const;
		double Convert_ADCToTime(const DRFDigiTime* locRFDigiTime) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double dRFBunchPeriod;
		double dTScale_FADC250;
		double dTScale_CAEN;
		map<DetectorSystem_t, double> dTimeOffsetMap;
		map<DetectorSystem_t, double> dTimeOffsetMap_TDCs;
};

#endif // _DRFTime_factory_
