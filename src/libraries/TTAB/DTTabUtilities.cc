// $Id$
//
//    File: DTTabUtilities.cc
// Created: Fri Apr  3 09:41:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#include "DTTabUtilities.h"

DTTabUtilities::DTTabUtilities(JEventLoop* locEventLoop)
{
	//import run-dependent info here (either from CCDB or event stream)

	// Get DF1TDCConfig's, put into map (will use to convert F1TDC digi time to ns)
		// Should be constant for a given run
	vector<const DF1TDCConfig*> locF1TDCConfigs;
	locEventLoop->Get(locF1TDCConfigs);
	for(size_t loc_i = 0; loc_i < locF1TDCConfigs.size(); ++loc_i)
		dF1TDCConfigMap[locF1TDCConfigs[loc_i]->rocid] = locF1TDCConfigs[loc_i];

	//Early Commissioning Data: Code & CCDB constants
	//BCAL, RF: none
	//SC, TAGH, TAGM: Uses F1TDC rollover
	//FDC: is bad (tdc_scale fixed at 0.115, no ref time or tframe). Use F1TDC rollover x 2
	//PSC: Uses F1TDC rollover; if run = 0 change run = 2012 (same as old hardcoded)

	// F1TDC tframe(ns) and rollover count
	map<string, int> tdc_parms;
	if(locEventLoop->GetCalib("/F1TDC/rollover", tdc_parms))
		jout << "Error loading /F1TDC/rollover !" << endl;

	map<string, int>::const_iterator locMapIterator = tdc_parms.find("tframe");
	dRolloverTimeWindowLength = (locMapIterator != tdc_parms.end()) ? double(tdc_parms["tframe"]) : std::numeric_limits<double>::quiet_NaN();

	locMapIterator = tdc_parms.find("count");
	dNumTDCTicksInRolloverTimeWindow = (locMapIterator != tdc_parms.end()) ? double(tdc_parms["count"]) : std::numeric_limits<double>::quiet_NaN();

	if(locEventLoop->GetJEvent().GetRunNumber() == 0) //PSC data with bad run number. Use hard-coded values from run 2012
	{
		dRolloverTimeWindowLength = 3744;
		dNumTDCTicksInRolloverTimeWindow = 64466;
	}
}

double DTTabUtilities::Convert_DigiTimeToNs(const JObject* locTDCDigiHit) const
{
	//See if the input object is an DF1TDCHit. If so, convert it
	const DF1TDCHit* locF1TDCHit = dynamic_cast<const DF1TDCHit*>(locTDCDigiHit);
	if(locF1TDCHit != NULL) //it's an F1TDCHit
		return Convert_DigiTimeToNs(locF1TDCHit);

	//Get the DF1TDCHit associated object
	vector<const DF1TDCHit*> locF1TDCHits;
	locTDCDigiHit->Get(locF1TDCHits);
	if(locF1TDCHits.empty())
	{
		cout << "ERROR: INCORRECT INPUT OBJECT TO DTTabUtilities::Convert_DigiTimeToNs(). RETURNING NaN." << endl;
		return std::numeric_limits<double>::quiet_NaN();
	}

	//Convert
	return Convert_DigiTimeToNs(locF1TDCHits[0]);
}

double DTTabUtilities::Convert_DigiTimeToNs(const DF1TDCHit* locF1TDCHit) const
{
    uint32_t locROCID = locF1TDCHit->rocid;

    // Get iterators for the DCODAROCInfo and DF1TDCConfig for this ROC
    map<uint32_t, const DCODAROCInfo*>::const_iterator locROCInfoIterator = dCODAROCInfoMap.find(locROCID);
    map<uint32_t, const DF1TDCConfig*>::const_iterator locF1TDCConfigIterator = dF1TDCConfigMap.find(locROCID);

    // Get DCODAROCInfo and DF1TDCConfig for this ROC
    const DF1TDCConfig* locF1TDCConfig = (locF1TDCConfigIterator != dF1TDCConfigMap.end()) ? locF1TDCConfigIterator->second : NULL;
    const DCODAROCInfo* locCODAROCInfo = (locROCInfoIterator != dCODAROCInfoMap.end()) ? locROCInfoIterator->second : NULL;

    if(locCODAROCInfo == NULL) //e.g. MC
    	return Convert_DigiTimeToNs_TriggerReferenceSignal(locF1TDCHit);

    if(locF1TDCConfig == NULL) //e.g. Early Fall 2014 Commissioning Data (use CCDB constants)
    	return Convert_DigiTimeToNs_GlobalSystemClock_CCDB(locF1TDCHit, locCODAROCInfo);

    // Have all objects needed, call the main function
	return Convert_DigiTimeToNs_GlobalSystemClock_ConfigInfo(locF1TDCHit, locCODAROCInfo, locF1TDCConfig);
}

double DTTabUtilities::Convert_DigiTimeToNs_GlobalSystemClock_ConfigInfo(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo, const DF1TDCConfig* locF1TDCConfig) const
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
	// DCODAROCInfo::timestamp is the number of clock ticks of the global system clock since it was reset (at the beginning of the run)
    uint64_t locReferenceClockTime = (uint64_t(4))*locCODAROCInfo->timestamp; // one clock tick = 4 ns
	uint64_t locNumRollovers = locReferenceClockTime/locRolloverTimeWindowLength;
	uint64_t locReferenceTimeThisWindow = locReferenceClockTime - locNumRollovers*locRolloverTimeWindowLength;

	//compute and return the time difference
	double locDeltaT = locTDCToNsScaleFactor*double(locF1TDCHit->time) - double(locReferenceTimeThisWindow); // in ns
	if(locDeltaT < 0.0) // Take care of rollover
		locDeltaT += double(locRolloverTimeWindowLength); //the time rolled over between the hit and reference times

	return locDeltaT;
}

double DTTabUtilities::Convert_DigiTimeToNs_GlobalSystemClock_CCDB(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo) const
{
	//compare the digi time to the global system clock available on the crate, then convert to ns
		//GlueX Doc 2686 Appendix A, using the global system clock, but the CCDB for the constants

    //GlueX DOC 747
    //48-channel-readout = F1TDCV3 (FDC) = low-resolution readout (115 ps vs 57 ps) //fewer (factor 2) TDC-ticks per ns
	bool locIsLowResolutionReadout = (locF1TDCHit->modtype == DModuleType::F1TDC48);

	//Calculate TDC -> ns Scale Factor
	double locTDCToNsScaleFactor = dRolloverTimeWindowLength/double(dNumTDCTicksInRolloverTimeWindow);
	if(locIsLowResolutionReadout) //should use HIGHRES bit from the data-stream, but it's not available
		locTDCToNsScaleFactor *= 2.0; //fewer TDC-ticks per ns

	// Transition the reference time into this F1TDC time window (get remainder)
	// DCODAROCInfo::timestamp is the number of clock ticks of the global system clock since it was reset (at the beginning of the run)
    uint64_t locReferenceClockTime = (uint64_t(4))*locCODAROCInfo->timestamp; // one clock tick = 4 ns
	uint64_t locNumRollovers = locReferenceClockTime/dRolloverTimeWindowLength;
	uint64_t locReferenceTimeThisWindow = locReferenceClockTime - locNumRollovers*dRolloverTimeWindowLength;

	//compute and return the time difference
	double locDeltaT = locTDCToNsScaleFactor*double(locF1TDCHit->time) - double(locReferenceTimeThisWindow); // in ns
	if(locDeltaT < 0.0) // Take care of rollover
		locDeltaT += double(dRolloverTimeWindowLength); //the time rolled over between the hit and reference times

	return locDeltaT;
}

double DTTabUtilities::Convert_DigiTimeToNs_TriggerReferenceSignal(const DF1TDCHit* locF1TDCHit) const
{
	//compare the digi time to the trigger reference signal, then convert to ns
		//GlueX Doc 2686 Appendix A, using the trigger reference signal

    //GlueX DOC 747
    //48-channel-readout = F1TDCV3 (FDC) = low-resolution readout (115 ps vs 57 ps) //fewer (factor 2) TDC-ticks per ns
	bool locIsLowResolutionReadout = (locF1TDCHit->modtype == DModuleType::F1TDC48);

	//Calculate TDC -> ns Scale Factor
	double locTDCToNsScaleFactor = dRolloverTimeWindowLength/double(dNumTDCTicksInRolloverTimeWindow);
	if(locIsLowResolutionReadout) //should use HIGHRES bit from the data-stream, but it's not available
		locTDCToNsScaleFactor *= 2.0; //fewer TDC-ticks per ns

	int64_t locDeltaTDC = int64_t(uint64_t(locF1TDCHit->time) - dTriggerReferenceSignal);

	// Take care of rollover
	if(locDeltaTDC < 0)
		locDeltaTDC += int64_t(dNumTDCTicksInRolloverTimeWindow);

	return locTDCToNsScaleFactor * double(locDeltaTDC);
}
