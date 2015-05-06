// $Id$
//
//    File: DTTabUtilities.cc
// Created: Fri Apr  3 09:41:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#include "DTTabUtilities.h"

DTTabUtilities::DTTabUtilities(void)
{
	dTScale_CAEN = 0.0234375; // ~ 23.4375 ps/count (TOF)
}

double DTTabUtilities::Convert_DigiTimeToNs_F1TDC(const JObject* locTDCDigiHit) const
{
	//See if the input object is an DF1TDCHit. If so, convert it
	const DF1TDCHit* locF1TDCHit = dynamic_cast<const DF1TDCHit*>(locTDCDigiHit);
	if(locF1TDCHit != NULL) //it's an F1TDCHit
		return Convert_DigiTimeToNs_F1TDC(locF1TDCHit);

	//Get the DF1TDCHit associated object
	vector<const DF1TDCHit*> locF1TDCHits;
	locTDCDigiHit->Get(locF1TDCHits);
	if(locF1TDCHits.empty())
	{
		cout << "ERROR: INCORRECT INPUT OBJECT TO DTTabUtilities::Convert_DigiTimeToNs_F1TDC(). RETURNING NaN." << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	//Convert
	return Convert_DigiTimeToNs_F1TDC(locF1TDCHits[0]);
}

double DTTabUtilities::Convert_DigiTimeToNs_F1TDC(const DF1TDCHit* locF1TDCHit) const
{
    uint32_t locROCID = locF1TDCHit->rocid;

    // Get DF1TDCConfig for this ROC
    vector<const DF1TDCConfig*> locF1TDCConfigs;
    locF1TDCHit->Get(locF1TDCConfigs);

    // Get DCODAROCInfo for this ROC
    map<uint32_t, const DCODAROCInfo*>::const_iterator locROCInfoIterator = dCODAROCInfoMap.find(locROCID);
    const DCODAROCInfo* locCODAROCInfo = (locROCInfoIterator != dCODAROCInfoMap.end()) ? locROCInfoIterator->second : NULL;

    if(locCODAROCInfo == NULL) //e.g. MC
    	return Convert_DigiTimeToNs_F1TDC_TriggerReferenceSignal(locF1TDCHit);

    if(locF1TDCConfigs.empty() || dHasBadOrNoF1TDCConfigInfoFlag) //e.g. Early Fall 2014 Commissioning Data (use CCDB constants)
    	return Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_CCDB(locF1TDCHit, locCODAROCInfo);

    // Have all objects needed, call the main function
	return Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_ConfigInfo(locF1TDCHit, locCODAROCInfo, locF1TDCConfigs[0]);
}

double DTTabUtilities::Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_ConfigInfo(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo, const DF1TDCConfig* locF1TDCConfig) const
{
	//compare the digi time to the global system clock available on the crate, then convert to ns
		//GlueX Doc 2686 Appendix A, using the global system clock

    //GlueX DOC 747
    //48-channel-readout = F1TDCV3 (FDC) = low-resolution readout (115 ps vs 57 ps) //fewer (factor 2) TDC-ticks per ns
	bool locIsLowResolutionReadout = (locF1TDCHit->modtype == DModuleType::F1TDC48);

	//Calculate TDC -> ns Scale Factor
	//32.0 is TREF and should be extracted from the data-stream but it's not available
	double locTDCToNsScaleFactor = (32.0/152.0) * double(locF1TDCConfig->REFCLKDIV) / (double(locF1TDCConfig->HSDIV)); // 32 ns / #-tdc-ticks-in-32ns
	if(locIsLowResolutionReadout) //should use HIGHRES bit from the data-stream, but it's not available
		locTDCToNsScaleFactor *= 2.0; //fewer TDC-ticks per ns

	//Calculate rollover window size for the F1TDCs
	uint64_t locRolloverTimeWindowLength = (uint64_t(32))*(uint64_t(locF1TDCConfig->REFCNT + 2));

	// Transition the reference time into this F1TDC time window (get remainder)
	uint64_t locReferenceTimeThisWindow = Calc_ROCRefTimeThisWindow(locCODAROCInfo, locRolloverTimeWindowLength);

	//compute and return the time difference
	double locDeltaT = locTDCToNsScaleFactor*double(locF1TDCHit->time) - double(locReferenceTimeThisWindow); // in ns
	if(locDeltaT < 0.0) // Take care of rollover
		locDeltaT += double(locRolloverTimeWindowLength); //the time rolled over between the hit and reference times

	return locDeltaT;
}

double DTTabUtilities::Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_CCDB(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo) const
{
	//compare the digi time to the global system clock available on the crate, then convert to ns
		//GlueX Doc 2686 Appendix A, using the global system clock, but the CCDB for the constants

    //GlueX DOC 747
    //48-channel-readout = F1TDCV3 (FDC) = low-resolution readout (115 ps vs 57 ps) //fewer (factor 2) TDC-ticks per ns
	bool locIsLowResolutionReadout = (locF1TDCHit->modtype == DModuleType::F1TDC48);

	//Calculate TDC -> ns Scale Factor
	double locTDCToNsScaleFactor = Calc_TDCToNsScaleFactor_CCDB(locIsLowResolutionReadout);

	// Transition the reference time into this F1TDC time window (get remainder)
	uint64_t locReferenceTimeThisWindow = Calc_ROCRefTimeThisWindow(locCODAROCInfo, dRolloverTimeWindowLength);

	//compute and return the time difference
	double locDeltaT = locTDCToNsScaleFactor*double(locF1TDCHit->time) - double(locReferenceTimeThisWindow); // in ns
	if(locDeltaT < 0.0) // Take care of rollover
		locDeltaT += double(dRolloverTimeWindowLength); //the time rolled over between the hit and reference times

	return locDeltaT;
}

double DTTabUtilities::Convert_DigiTimeToNs_F1TDC_TriggerReferenceSignal(const DF1TDCHit* locF1TDCHit) const
{
	//compare the digi time to the trigger reference signal, then convert to ns
		//GlueX Doc 2686 Appendix A, using the trigger reference signal

    //GlueX DOC 747
    //48-channel-readout = F1TDCV3 (FDC) = low-resolution readout (115 ps vs 57 ps) //fewer (factor 2) TDC-ticks per ns
	bool locIsLowResolutionReadout = (locF1TDCHit->modtype == DModuleType::F1TDC48);

	//Calculate TDC -> ns Scale Factor
	double locTDCToNsScaleFactor = Calc_TDCToNsScaleFactor_CCDB(locIsLowResolutionReadout);

	//Calc delta-t //potentially different scale factors between trigger reference & hit
	double locDeltaT = locTDCToNsScaleFactor*double(locF1TDCHit->time) - Convert_TriggerReferenceSignal();

	// Take care of rollover
	if(locDeltaT < 0)
		locDeltaT += double(dRolloverTimeWindowLength);

	return locDeltaT;
}

double DTTabUtilities::Calc_TDCToNsScaleFactor_CCDB(bool locIsLowResolutionReadout) const
{
	double locTDCToNsScaleFactor = double(dRolloverTimeWindowLength)/double(dNumTDCTicksInRolloverTimeWindow);
	if(locIsLowResolutionReadout) //should use HIGHRES bit from the data-stream, but it's not available
		locTDCToNsScaleFactor *= 2.0; //fewer TDC-ticks per ns
	return locTDCToNsScaleFactor;
}

uint64_t DTTabUtilities::Calc_ROCRefTimeThisWindow(const DCODAROCInfo* locCODAROCInfo, uint64_t locRolloverTimeWindowLength) const
{
	// Transition the reference time into this F1TDC time window (get remainder)
	// DCODAROCInfo::timestamp is the number of clock ticks of the global system clock since it was reset (at the beginning of the run)
    uint64_t locReferenceClockTime = (uint64_t(4))*locCODAROCInfo->timestamp; // one clock tick = 4 ns
	uint64_t locNumRollovers = locReferenceClockTime/locRolloverTimeWindowLength;
	return (locReferenceClockTime - locNumRollovers*locRolloverTimeWindowLength);
}

double DTTabUtilities::Convert_TriggerReferenceSignal(void) const
{
	double locTDCToNsScaleFactor = Calc_TDCToNsScaleFactor_CCDB(dTriggerReferenceSignalIsLowResTDC);
	return locTDCToNsScaleFactor * double(dTriggerReferenceSignal);
}

double DTTabUtilities::Convert_DigiTimeToNs_CAEN1290TDC(const JObject* locTDCDigiHit) const
{
	//See if the input object is an DCAEN1290TDCHit. If so, convert it
	const DCAEN1290TDCHit* locCAEN1290TDCHit = dynamic_cast<const DCAEN1290TDCHit*>(locTDCDigiHit);
	if(locCAEN1290TDCHit != NULL) //it's an DCAEN1290TDCHit
		return Convert_DigiTimeToNs_CAEN1290TDC(locCAEN1290TDCHit);

	//Get the DF1TDCHit associated object
	vector<const DCAEN1290TDCHit*> locCAEN1290TDCHits;
	locTDCDigiHit->Get(locCAEN1290TDCHits);
	if(locCAEN1290TDCHits.empty())
	{
		cout << "ERROR: INCORRECT INPUT OBJECT TO DTTabUtilities::Convert_DigiTimeToNs_CAEN1290TDC(). RETURNING NaN." << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	//Convert
	return Convert_DigiTimeToNs_CAEN1290TDC(locCAEN1290TDCHits[0]);
}

double DTTabUtilities::Convert_DigiTimeToNs_CAEN1290TDC(const DCAEN1290TDCHit* locCAEN1290TDCHit) const
{
	//compare the digi time to the trigger reference signal, then convert to ns
		//GlueX Doc 2686
    uint32_t locROCID = locCAEN1290TDCHit->rocid;

    // Get DCODAROCInfo for this ROC
    map<uint32_t, const DCODAROCInfo*>::const_iterator locROCInfoIterator = dCODAROCInfoMap.find(locROCID);
    if(locROCInfoIterator == dCODAROCInfoMap.end())
    	return std::numeric_limits<double>::quiet_NaN();
    const DCODAROCInfo* locCODAROCInfo = locROCInfoIterator->second;

    int locSystemClockBinRemainder = locCODAROCInfo->timestamp % 6;

	// The number of TI-counter (4 ns) blocks to shift the CAEN time to line-up with the TI time
    int locNum4NsBlocksToShift = dCAENTIPhaseDifference - locSystemClockBinRemainder;
    if(locNum4NsBlocksToShift < 0)
    	locNum4NsBlocksToShift += 6;

    // TDC bins are 25 ps wide, so each 4ns-TI-block is 160 bins ("TDC-ticks")
    int locTDCShift = 160 * locNum4NsBlocksToShift;

    return dTScale_CAEN*double(locCAEN1290TDCHit->time + locTDCShift);
}
