// $Id$
//
//    File: DRFTime_factory_FDC.h
// Created: Mon Mar 30 10:51:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DRFTime_factory_FDC_
#define _DRFTime_factory_FDC_

#include <limits>
#include <iostream>

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DRFTime.h"
#include "DRFDigiTime.h"
#include "DRFTDCDigiTime.h"

using namespace std;
using namespace jana;

class DRFTime_factory_FDC : public jana::JFactory<DRFTime>
{
	public:
		DRFTime_factory_FDC(){};
		~DRFTime_factory_FDC(){};
		const char* Tag(void){return "FDC";}

		double Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo) const;
		double Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo, double locPeriod) const;

		double Convert_TDCToTime(const DRFTDCDigiTime* locRFTDCDigiTime, const DTTabUtilities* locTTabUtilities) const;

		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
	
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double Calc_WeightedAverageRFTime(vector<double>& locRFTimes, double& locRFTimeVariance) const;

		DetectorSystem_t dSourceSystem;
		double dBeamBunchPeriod;

		double dTimeOffset;
		double dTimeOffsetVariance;
		double dTimeResolutionSq;
};

#endif // _DRFTime_factory_FDC_

