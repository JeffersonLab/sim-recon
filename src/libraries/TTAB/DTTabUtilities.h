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

using namespace std;

// This object is MODIFIED every event, so make sure to get it anew for each event!
	//DO NOT GRAB in your factory's brun() method.
class DTTabUtilities : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DTTabUtilities);
		
		DTTabUtilities(JEventLoop* locEventLoop);
		double Convert_DigiTimeToNs(const JObject* locTDCDigiHit) const;

		void Set_TriggerReferenceSignal(uint64_t locTriggerReferenceSignal){dTriggerReferenceSignal = locTriggerReferenceSignal;}
		void Set_CODAROCInfoMap(const map<uint32_t, const DCODAROCInfo*>& locCODAROCInfoMap){dCODAROCInfoMap = locCODAROCInfoMap;}

	private:
		//New System
		map<uint32_t, const DCODAROCInfo*> dCODAROCInfoMap; //key is rocid
		map<uint32_t, const DF1TDCConfig*> dF1TDCConfigMap; //key is rocid

		//Old System ONLY //Early Fall 2014 Commissioning data ONLY
		uint64_t dTriggerReferenceSignal;
		double dRolloverTimeWindowLength; //"T" or "T_{frame}"
		uint64_t dNumTDCTicksInRolloverTimeWindow; //"N" or "N_{frame}"

		double Convert_DigiTimeToNs(const DF1TDCHit* locF1TDCHit) const;
		double Convert_DigiTimeToNs_GlobalSystemClock_ConfigInfo(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo, const DF1TDCConfig* locF1TDCConfig) const;
		double Convert_DigiTimeToNs_GlobalSystemClock_CCDB(const DF1TDCHit* locF1TDCHit, const DCODAROCInfo* locCODAROCInfo) const;
		double Convert_DigiTimeToNs_TriggerReferenceSignal(const DF1TDCHit* locF1TDCHit) const;
};

#endif // _DTTabUtilities_
