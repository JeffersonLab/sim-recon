// $Id$
//
//    File: DTTabUtilities.h
// Created: Fri Apr  3 09:41:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DTTabUtilities_
#define _DTTabUtilities_

#include <limits>
#include <map>
#include <string>
#include <vector>

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include <DAQ/DCODAROCInfo.h>
#include <DAQ/DF1TDCConfig.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DCAEN1290TDCHit.h>

using namespace std;

// This object is MODIFIED every event, so make sure to get it anew for each event!
	//DO NOT GRAB in your factory's brun() method.
class DTTabUtilities : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DTTabUtilities);
		
		DTTabUtilities(void);

		double Convert_DigiTimeToNs_F1TDC(const JObject* locTDCDigiHit) const;
		double Convert_DigiTimeToNs_F1TDC(const DF1TDCHit* locF1TDCHit) const;
		double Convert_DigiTimeToNs_CAEN1290TDC(const JObject* locTDCDigiHit) const;
		double Convert_DigiTimeToNs_CAEN1290TDC(const DCAEN1290TDCHit* locCAEN1290TDCHit) const;

		//F1TDCs: New System
		bool dHasBadOrNoF1TDCConfigInfoFlag;
		map<uint32_t, const DCODAROCInfo*> dCODAROCInfoMap; //key is rocid

		//F1TDCs: Old System ONLY //Early Fall 2014 Commissioning data ONLY
		uint64_t dTriggerReferenceSignal;
		bool dTriggerReferenceSignalIsLowResTDC;
		uint64_t dRolloverTimeWindowLength; //"T" or "T_{frame}"
		uint64_t dNumTDCTicksInRolloverTimeWindow; //"N" or "N_{frame}"

		//CAEN1290s:
		int dCAENTIPhaseDifference; //0 -> 5
		double dTScale_CAEN;

		uint64_t Calc_ROCRefTimeThisWindow(const DCODAROCInfo* locCODAROCInfo, uint64_t locRolloverTimeWindowLength) const;
		double Calc_TDCToNsScaleFactor_CCDB(bool locIsLowResolutionReadout) const;
		double Convert_TriggerReferenceSignal(void) const;

        //FADC250s: convenience functions for v2 firmware (fall 2016 -> ?)
        bool CheckFADC250_PedestalOK(uint32_t QF) const;
        bool CheckFADC250_NoErrors(uint32_t QF) const;

	private:

		double Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_ConfigInfo(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo, const DF1TDCConfig* locF1TDCConfig) const;
		double Convert_DigiTimeToNs_F1TDC_GlobalSystemClock_CCDB(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo) const;
		double Convert_DigiTimeToNs_F1TDC_TriggerReferenceSignal(const DF1TDCHit* locF1TDCHit) const;
};

#endif // _DTTabUtilities_
