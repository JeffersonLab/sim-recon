/*
 * RICH_init.cc
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DRichHit.h"

typedef JFactory<DRichHit> DRichHit_factory;

jerror_t RICH_init(JEventLoop *loop) {
	/// Create and register RICH data factories
	loop->AddFactory(new DRichHit_factory());

	return NOERROR;
}



