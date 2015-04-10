// $Id$
//
//    File: DRFTime_factory.cc
// Created: Mon Mar 30 10:51:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#include "DRFTime_factory.h"

//------------------
// init
//------------------
jerror_t DRFTime_factory::init(void)
{
	// These are (apparently) fixed for all time
	dTScale_FADC250 = 0.0625; // ns per count
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DRFTime_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	map<string, double> locCCDBMap;
	map<string, double>::iterator locIterator;
	DetectorSystem_t locSystem;
	bool locIsTDCFlag;

	//RF Period
	vector<double> locRFPeriodVector;
	locEventLoop->GetCalib("PHOTON_BEAM/RF/rf_period", locRFPeriodVector);
	dRFBunchPeriod = locRFPeriodVector[0];

	//Time Offsets
	if(locEventLoop->GetCalib("/PHOTON_BEAM/RF/time_offset", locCCDBMap))
		jout << "Error loading /PHOTON_BEAM/RF/time_offset !" << endl;
	for(locIterator = locCCDBMap.begin(); locIterator != locCCDBMap.end(); ++locIterator)
	{
		Extract_DetectorSystemAndType(locIterator->first, locSystem, locIsTDCFlag);
		map<DetectorSystem_t, double>& locFactoryMap = locIsTDCFlag ? dTimeOffsetMap_TDCs : dTimeOffsetMap_ADCs;
		locFactoryMap[locSystem] = locIterator->second;
	}

	//Time Offset Variances
	if(locEventLoop->GetCalib("/PHOTON_BEAM/RF/time_offset_var", locCCDBMap))
		jout << "Error loading /PHOTON_BEAM/RF/time_offset_var !" << endl;
	for(locIterator = locCCDBMap.begin(); locIterator != locCCDBMap.end(); ++locIterator)
	{
		Extract_DetectorSystemAndType(locIterator->first, locSystem, locIsTDCFlag);
		map<DetectorSystem_t, double>& locFactoryMap = locIsTDCFlag ? dTimeOffsetVarianceMap_TDCs : dTimeOffsetVarianceMap_ADCs;
		locFactoryMap[locSystem] = locIterator->second;
	}

	//Time Resolution Squared
	if(locEventLoop->GetCalib("/PHOTON_BEAM/RF/time_resolution_sq", locCCDBMap))
		jout << "Error loading /PHOTON_BEAM/RF/time_resolution_sq !" << endl;
	for(locIterator = locCCDBMap.begin(); locIterator != locCCDBMap.end(); ++locIterator)
	{
		Extract_DetectorSystemAndType(locIterator->first, locSystem, locIsTDCFlag);
		map<DetectorSystem_t, double>& locFactoryMap = locIsTDCFlag ? dTimeResolutionSqMap_TDCs : dTimeResolutionSqMap_ADCs;
		locFactoryMap[locSystem] = locIterator->second;
	}

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

	//Get All RF Times
	map<pair<DetectorSystem_t, bool>, vector<double> > locRFTimesMap; //bool = isTDC

	//F1TDC's
	for(size_t loc_i = 0; loc_i < locRFTDCDigiTimes.size(); ++loc_i)
	{
		const DRFTDCDigiTime* locRFTDCDigiTime = locRFTDCDigiTimes[loc_i];
		DetectorSystem_t locSystem = locRFTDCDigiTime->dSystem;
		if(dTimeOffsetMap_TDCs.find(locSystem) == dTimeOffsetMap_TDCs.end())
			continue; //system not recognized
		double locRFTime = Convert_TDCToTime(locRFTDCDigiTime, locTTabUtilities) - dTimeOffsetMap_TDCs[locSystem];
		locRFTimesMap[pair<DetectorSystem_t, bool>(locSystem, true)].push_back(locRFTime);
	}

	//FADC250s
	for(size_t loc_i = 0; loc_i < locRFDigiTimes.size(); ++loc_i)
	{
		const DRFDigiTime* locRFDigiTime = locRFDigiTimes[loc_i];
		DetectorSystem_t locSystem = locRFDigiTime->dSystem;
		if(dTimeOffsetMap_ADCs.find(locSystem) == dTimeOffsetMap_ADCs.end())
			continue; //system not recognized
		double locRFTime = Convert_ADCToTime(locRFDigiTime) - dTimeOffsetMap_ADCs[locSystem];
		locRFTimesMap[pair<DetectorSystem_t, bool>(locSystem, false)].push_back(locRFTime);
	}

	if(locRFTimesMap.empty())
		return NOERROR; //No RF signals, will try to emulate RF bunch time from timing from other systems

	//take weighted average of times
		//avg = Sum(w_i * x_i) / Sum (w_i)
			//w_i = 1/varx_i
		//avg = Sum(x_i / varx_i) / Sum (1 / varx_i)
			//var_avg = 1 / Sum (1 / varx_i)
		//avg = Sum(x_i / varx_i) * var_avg

	//will line up the first time near zero, then line up subsequent times to the first time
		//cannot line up all times to zero: if first time is near +/- rf_period/2 (e.g. +/- 1.002 ns),
		//may end up with some times near 1.002 and some near -1.002!

	map<pair<DetectorSystem_t, bool>, vector<double> >::iterator locTimeIterator = locRFTimesMap.begin();
	double locFirstTime = (locTimeIterator->second)[0];
	locFirstTime = Step_TimeToNearInputTime(locFirstTime, 0.0);

	double locSumOneOverTimeVariance = 0.0;
	double locSumTimeOverTimeVariance = 0.0;
	for(; locTimeIterator != locRFTimesMap.end(); ++locTimeIterator)
	{
		pair<DetectorSystem_t, bool> locSystemPair = locTimeIterator->first;
		vector<double>& locRFTimes = locTimeIterator->second;
		DetectorSystem_t locSystem = locSystemPair.first;
		bool locIsTDCFlag = locSystemPair.second;

		double locSingleTimeVariance = 0.0;
		if(locIsTDCFlag)
			locSingleTimeVariance = dTimeResolutionSqMap_TDCs[locSystem] + dTimeOffsetVarianceMap_TDCs[locSystem];
		else
			locSingleTimeVariance = dTimeResolutionSqMap_ADCs[locSystem] + dTimeOffsetVarianceMap_ADCs[locSystem];
		locSumOneOverTimeVariance += double(locRFTimes.size()) / locSingleTimeVariance;

		for(size_t loc_i = 0; loc_i < locRFTimes.size(); ++loc_i)
		{
			double locShiftedRFTime = Step_TimeToNearInputTime(locRFTimes[loc_i], locFirstTime);
			locSumTimeOverTimeVariance += locShiftedRFTime / locSingleTimeVariance;
		}
	}
	double locRFTimeVariance = 1.0 / locSumOneOverTimeVariance;
	double locAverageRFTime = locSumTimeOverTimeVariance * locRFTimeVariance;

	DRFTime* locRFTime = new DRFTime();
	locRFTime->dTime = locAverageRFTime; //This time is defined at the center of the target (CCDB offsets center it)
	locRFTime->dTimeVariance = locRFTimeVariance;
	_data.push_back(locRFTime);

	return NOERROR;
}

void DRFTime_factory::Extract_DetectorSystemAndType(string locKeyName, DetectorSystem_t& locSystem, bool& locIsTDCFlag) const
{
	//Assumes Input is in form "System_Type" (e.g. "TOF_TDC", "TAGH_ADC")
	size_t locUnderscoreIndex = locKeyName.find_first_of("_");
	string locSystemName = locKeyName.substr(0, locUnderscoreIndex);
	locSystem = NameToSystem(locSystemName.c_str());
	string locTypeName = locKeyName.substr(locUnderscoreIndex + 1);
	locIsTDCFlag = (locTypeName == "TDC");
}

double DRFTime_factory::Step_TimeToNearInputTime(double locTimeToStep, double locTimeToStepTo) const
{
	double locDeltaT = locTimeToStepTo - locTimeToStep;
	int locNumRFBucketsToShift = (locDeltaT > 0.0) ? int(locDeltaT/dRFBunchPeriod + 0.5) : int(locDeltaT/dRFBunchPeriod - 0.5);
	return (locTimeToStep + dRFBunchPeriod*double(locNumRFBucketsToShift));
}

double DRFTime_factory::Convert_TDCToTime(const DRFTDCDigiTime* locRFTDCDigiTime, const DTTabUtilities* locTTabUtilities) const
{
	if(locRFTDCDigiTime->dIsCAENTDCFlag)
		return locTTabUtilities->Convert_DigiTimeToNs_CAEN1290TDC(locRFTDCDigiTime);
	else
		return locTTabUtilities->Convert_DigiTimeToNs_F1TDC(locRFTDCDigiTime);
}

double DRFTime_factory::Convert_ADCToTime(const DRFDigiTime* locRFDigiTime) const
{
	return dTScale_FADC250*locRFDigiTime->time;
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
