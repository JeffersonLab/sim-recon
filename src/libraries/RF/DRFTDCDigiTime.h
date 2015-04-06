// $Id$
//
//    File: DRFTDCDigiTime.h
// Created: Mon Mar 30 10:51:39 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DRFTDCDigiTime_
#define _DRFTDCDigiTime_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include "DAQ/DModuleType.h"
#include "GlueX.h"

using namespace std;

class DRFTDCDigiTime : public jana::JObject
{
	public:
		JOBJECT_PUBLIC(DRFTDCDigiTime);

		DetectorSystem_t dSystem;
		bool dIsCAENTDCFlag;
		
		uint32_t time;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "System", "%s", SystemName(dSystem));
			AddString(items, "time", "%u", time);
		}
};

#endif // _DRFTDCDigiTime_

