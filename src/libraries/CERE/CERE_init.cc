/*
 * CERE_init.cc
 *
 *  Created on: Oct 3, 2012
 *      Author: yqiang
 *
 *  Modified on:
 *  	Oct 9 2012, Yi Qiang, add general Cerenkov hit
 */

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DCereHit.h"

jerror_t CERE_init(JEventLoop *loop) {
	/// Create and register Cherenkov data factories
	loop->AddFactory(new JFactory<DCereHit>());

	return NOERROR;
}

