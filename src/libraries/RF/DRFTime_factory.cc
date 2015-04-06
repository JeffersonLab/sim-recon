// $Id$
//
//    File: DRFTime_factory.cc
// Created: Mon Mar 30 10:51:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DRFTime_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DRFTime_factory::init(void)
{
	// These are (apparently) fixed for all time
	dTScale_FADC250 = 0.0625; // ns per count
	dTScale_CAEN = 0.025; // 25 ps/count (TOF)
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DRFTime_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	dRFBunchPeriod = 2.004; //GET FROM CCDB

	map<DetectorSystem_t, double> dTimeOffsetMap; //GET FROM CCDB
	map<DetectorSystem_t, double> dTimeOffsetMap_TDCs; //GET FROM CCDB
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DRFTime_factory::evnt(JEventLoop *locEventLoop, int eventnumber)
{
	//The RF Time is defined at the center of the target
	//The time offset needed to line it up to the center of the target is absorbed into the calibration constants.
	vector<const DRFDigiTime*> locRFDigiTimes;
	locEventLoop->Get(locRFDigiTimes);

	vector<const DRFTDCDigiTime*> locRFTDCDigiTimes;
	locEventLoop->Get(locRFTDCDigiTimes);

	const DTTabUtilities* locTTabUtilities = NULL;
	locEventLoop->GetSingle(locTTabUtilities);

	//choose one RF time object (time) from each source
	map<DetectorSystem_t, double> locRFTimeBySystemMap;

	//F1TDC's
	for(size_t loc_i = 0; loc_i < locRFTDCDigiTimes.size(); ++loc_i)
	{
		const DRFTDCDigiTime* locRFTDCDigiTime = locRFTDCDigiTimes[loc_i];
		DetectorSystem_t locSystem = locRFTDCDigiTime->dSystem;
		if(locRFTimeBySystemMap.find(locSystem) != locRFTimeBySystemMap.end())
			continue; //time already found for this system
		if(dTimeOffsetMap_TDCs.find(locSystem) == dTimeOffsetMap_TDCs.end())
			continue; //system not recognized

		double locRFTime = -1.0*dTimeOffsetMap_TDCs[locSystem];
		if(locRFTDCDigiTime->dIsCAENTDCFlag)
			locRFTime += dTScale_CAEN*double(locRFTDCDigiTime->time);
		else
			locRFTime += locTTabUtilities->Convert_DigiTimeToNs(locRFTDCDigiTime);

		locRFTimeBySystemMap[locSystem] = locRFTime;
	}

	//FADC250's if no F1TDC's (have worse resolution)
	if(locRFTimeBySystemMap.empty())
	{
		for(size_t loc_i = 0; loc_i < locRFDigiTimes.size(); ++loc_i)
		{
			const DRFDigiTime* locRFDigiTime = locRFDigiTimes[loc_i];
			DetectorSystem_t locSystem = locRFDigiTime->dSystem;
			if(locRFTimeBySystemMap.find(locSystem) != locRFTimeBySystemMap.end())
				continue; //time already found for this system
			if(dTimeOffsetMap.find(locSystem) == dTimeOffsetMap.end())
				continue; //system not recognized
			double locRFTime = dTScale_FADC250*locRFDigiTime->time - dTimeOffsetMap[locSystem];
			locRFTimeBySystemMap[locSystem] = locRFTime;
		}
	}

	if(locRFTimeBySystemMap.empty())
		return NOERROR; //No RF signals, will try to emulate RF bunch time from timing from other systems

	//line up the times near zero, take the average
	map<DetectorSystem_t, double>::iterator locIterator = locRFTimeBySystemMap.begin();
	double locAverageRFTime = Step_TimeToNearInputTime(locIterator->second, 0.0);
	for(++locIterator; locIterator != locRFTimeBySystemMap.end(); ++locIterator)
		locAverageRFTime += Step_TimeToNearInputTime(locIterator->second, locAverageRFTime);
	locAverageRFTime /= double(locRFTimeBySystemMap.size());

	DRFTime* locRFTime = new DRFTime();
	locRFTime->dTime = locAverageRFTime; //This time is defined at the center of the target.
	locRFTime->dTimeVariance = 0.0; //SET ME!!!
	_data.push_back(locRFTime);

	return NOERROR;
}

double DRFTime_factory::Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo)
{
	while((locTimeToStep - locTimeToStepTo) > 0.5*dRFBunchPeriod)
		locTimeToStep -= dRFBunchPeriod;
	while((locTimeToStepTo - locTimeToStep) > 0.5*dRFBunchPeriod)
		locTimeToStep += dRFBunchPeriod;
	return locTimeToStep;
}

//------------------
// erun
//------------------
jerror_t DRFTime_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DRFTime_factory::fini(void)
{
	return NOERROR;
}

