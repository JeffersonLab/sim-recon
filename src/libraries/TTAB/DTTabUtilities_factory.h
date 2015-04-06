// $Id$
//
//    File: DTTabUtilities_factory.h
// Created: Fri Apr  3 09:41:29 EDT 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.12.2.el6.x86_64 x86_64)
//

#ifndef _DTTabUtilities_factory_
#define _DTTabUtilities_factory_

#include <JANA/JFactory.h>
#include "DTTabUtilities.h"

using namespace jana;

class DTTabUtilities_factory : public jana::JFactory<DTTabUtilities>
{
	public:
		DTTabUtilities_factory(){};
		virtual ~DTTabUtilities_factory(){};

	private:
		jerror_t brun(jana::JEventLoop *loop, int runnumber);
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);

		//New system
		map<uint32_t, const DF1TDCConfig*> dF1TDCConfigMap; //key is rocid

		//Old System ONLY //Early Fall 2014 Commissioning data ONLY
		uint64_t dRolloverTimeWindowLength; //"T" or "T_{frame}"
		uint64_t dNumTDCTicksInRolloverTimeWindow; //"N" or "N_{frame}"
};

#endif // _DTTabUtilities_factory_
