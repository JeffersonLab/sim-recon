// $Id$
//
//    File: DTTabUtilities_factory.cc
// Created: Mon Apr  6 09:41:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#include "DTTabUtilities_factory.h"

jerror_t DTTabUtilities_factory::brun(jana::JEventLoop* locEventLoop, int runnumber)
{
	// Get DF1TDCConfig's, put into map (will use to convert F1TDC digi time to ns)
		// Should be constant for a given run
	vector<const DF1TDCConfig*> locF1TDCConfigs;
	locEventLoop->Get(locF1TDCConfigs);

	dF1TDCConfigMap.clear();
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
	dRolloverTimeWindowLength = (locMapIterator != tdc_parms.end()) ? uint64_t(tdc_parms["tframe"]) : 0;

	locMapIterator = tdc_parms.find("count");
	dNumTDCTicksInRolloverTimeWindow = (locMapIterator != tdc_parms.end()) ? uint64_t(tdc_parms["count"]) : 0;

	if(locEventLoop->GetJEvent().GetRunNumber() == 0) //PSC data with bad run number. Use hard-coded values from run 2012
	{
		dRolloverTimeWindowLength = 3744;
		dNumTDCTicksInRolloverTimeWindow = 64466;
	}

	//CAEN1290/TI Phase Difference
	dCAENTIPhaseDifference = 1;
	map<string, double> tof_tdc_shift;
	if(!eventLoop->GetCalib("/TOF/tdc_shift", tof_tdc_shift))
		dCAENTIPhaseDifference = tof_tdc_shift["TOF_TDC_SHIFT"];
	dCAENTIPhaseDifference = 1;
	return NOERROR;
}

jerror_t DTTabUtilities_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	DTTabUtilities* locTTabUtilities = new DTTabUtilities();
	locTTabUtilities->dF1TDCConfigMap = dF1TDCConfigMap;
	locTTabUtilities->dRolloverTimeWindowLength = dRolloverTimeWindowLength;
	locTTabUtilities->dNumTDCTicksInRolloverTimeWindow = dNumTDCTicksInRolloverTimeWindow;
	locTTabUtilities->dIsFallCommissioningDataFlag = (locEventLoop->GetJEvent().GetRunNumber() <= 2700);
	locTTabUtilities->dCAENTIPhaseDifference = dCAENTIPhaseDifference;

	// Get DCODAROCInfo's, put into map
	vector<const DCODAROCInfo*> locCODAROCInfos;
	locEventLoop->Get(locCODAROCInfos);

	map<uint32_t, const DCODAROCInfo*> locCODAROCInfoMap;
	for(size_t loc_i = 0; loc_i < locCODAROCInfos.size(); ++loc_i)
		locCODAROCInfoMap[locCODAROCInfos[loc_i]->rocid] = locCODAROCInfos[loc_i];
	locTTabUtilities->dCODAROCInfoMap = locCODAROCInfoMap;

	//get the trigger reference signal ("Beni-cable")
		//hard-coded crate/slot/channel, but whatever. This isn't intended to be long-term-code anyway.

	vector<const DF1TDCHit*> locF1TDCHits;
	locEventLoop->Get(locF1TDCHits);

	bool locFoundFlag = false;
	for(size_t loc_i = 0; loc_i < locF1TDCHits.size(); ++loc_i)
	{
		if((locF1TDCHits[loc_i]->rocid != 51) || (locF1TDCHits[loc_i]->slot != 17) || (locF1TDCHits[loc_i]->channel != 8))
			continue;
		locTTabUtilities->dTriggerReferenceSignal = locF1TDCHits[loc_i]->time; //in TDC clicks
		locTTabUtilities->dTriggerReferenceSignalIsLowResTDC = (locF1TDCHits[loc_i]->modtype == DModuleType::F1TDC48);
		locFoundFlag = true;
		break;
	}
	if(!locFoundFlag)
		locTTabUtilities->dTriggerReferenceSignal = 0;

	_data.push_back(locTTabUtilities);
	return NOERROR;
}
