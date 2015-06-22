/*
 * DIRC_init.cc
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 *  Modified on: 0ct 7, 2013, yqiang, added DIRCTruthHit factory
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DDIRCHit.h"
#include "DDIRCTruthHit.h"

jerror_t DIRC_init(JEventLoop *loop) {
	/// Create and register DIRC data factories
	loop->AddFactory(new JFactory<DDIRCHit>());
	loop->AddFactory(new JFactory<DDIRCTruthHit>());

	return NOERROR;
}

