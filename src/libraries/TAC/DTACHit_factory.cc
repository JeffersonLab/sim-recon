/*
 * DTACHIT_factory.cc
 *
 *  Created on: Mar 24, 2017
 *      Author: hovanes
 */

#include "DTACHit_factory.h"

//------------------
// init
//------------------
jerror_t DTACHit_factory::init(void) {
	// Set default configuration parameters

	// Setting this flag makes it so that JANA does not delete the objects in _data.
	// This factory will manage this memory.
	SetFactoryFlag(NOT_OBJECT_OWNER);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTACHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber) {
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTACHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber) {
	// Clear _data vector
	Reset_Data();

	// Get TAC hits
	vector<const DTACHit*> hits;
	loop->Get(hits);

	for (auto& hit : hits) {
		if (hit->isFADCPresent()) {
			_data.push_back(const_cast<DTACHit*>(hit));
		}
	}
	return NOERROR;
}

//------------------
// Reset_Data()
//------------------
void DTACHit_factory::Reset_Data(void)
{
	// Clear _data vector
	_data.clear();
}

//------------------
// erun
//------------------
jerror_t DTACHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTACHit_factory::fini(void)
{
    return NOERROR;
}
