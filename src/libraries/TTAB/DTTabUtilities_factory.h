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
		jerror_t erun(void);
};

inline jerror_t DTTabUtilities_factory::brun(jana::JEventLoop *loop, int runnumber)
{
	// Create single DTTabUtilities object and mark the factory as
	// persistent so it doesn't get deleted every event.

	// However, it will get MODIFIED every event. Make sure to get it anew for each event!

	DTTabUtilities *locTTabUtilities = new DTTabUtilities(loop);
	SetFactoryFlag(PERSISTANT);
	ClearFactoryFlag(WRITE_TO_OUTPUT);
	_data.push_back(locTTabUtilities);

	return NOERROR;
}

inline jerror_t DTTabUtilities_factory::evnt(jana::JEventLoop *eventLoop, int eventnumber)
{
	// Get DCODAROCInfo's, put into map
	vector<const DCODAROCInfo*> locCODAROCInfos;
	eventLoop->Get(locCODAROCInfos);

	map<uint32_t, const DCODAROCInfo*> locCODAROCInfoMap;
	for(size_t loc_i = 0; loc_i < locCODAROCInfos.size(); ++loc_i)
		locCODAROCInfoMap[locCODAROCInfos[loc_i]->rocid] = locCODAROCInfos[loc_i];
	_data[0]->Set_CODAROCInfoMap(locCODAROCInfoMap);

	//get the trigger reference signal ("Beni-cable")
		//hard-coded crate/slot/channel, but whatever. This isn't intended to be long-term-code anyway.

	vector<const DF1TDCHit*> locF1TDCHits;
	eventLoop->Get(locF1TDCHits);

	bool locFoundFlag = false;
	for(size_t loc_i = 0; loc_i < locF1TDCHits.size(); ++loc_i)
	{
		if((locF1TDCHits[loc_i]->rocid != 51) || (locF1TDCHits[loc_i]->slot != 17) || (locF1TDCHits[loc_i]->channel != 8))
			continue;
		_data[0]->Set_TriggerReferenceSignal(locF1TDCHits[loc_i]->time); //in TDC clicks
		locFoundFlag = true;
		break;
	}
	if(!locFoundFlag)
		_data[0]->Set_TriggerReferenceSignal(std::numeric_limits<double>::quiet_NaN());

	return NOERROR;
}

inline jerror_t DTTabUtilities_factory::erun(void)
{
	delete _data[0];
	_data.clear();

	return NOERROR;
}

#endif // _DTTabUtilities_factory_
