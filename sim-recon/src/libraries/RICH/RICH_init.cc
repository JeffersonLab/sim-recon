/*
 * RICH_init.cc
 *
 *  Created on: Oct 11, 2012
 *      Author: yqiang
 *  Modified on: 0ct 7, 2013, yqiang, added RichTruthHit factory
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DRichHit.h"
#include "DRichTruthHit.h"

typedef JFactory<DRichHit> DRichHit_factory;
typedef JFactory<DRichTruthHit> DRichTruthHit_factory;

jerror_t RICH_init(JEventLoop *loop) {
	/// Create and register RICH data factories
	loop->AddFactory(new DRichHit_factory());
	loop->AddFactory(new DRichTruthHit_factory());

	return NOERROR;
}

