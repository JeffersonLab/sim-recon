/*
 * TAC_init.cc
 *
 *  Created on: Mar 28, 2017
 *      Author: hovanes
 */



#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTACDigiHit.h"
#include "DTACTDCDigiHit.h"
#include "DTACHit_factory.h"

jerror_t TAC_init(JEventLoop *loop)
{
	/// Create and register TAC data factories
	loop->AddFactory(new JFactory<DTACDigiHit>());
	loop->AddFactory(new JFactory<DTACTDCDigiHit>());
	loop->AddFactory(new DTACHit_factory());
	loop->AddFactory(new JFactory<DTACHit>("TRUTH"));

	return NOERROR;
}

