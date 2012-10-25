// $Id$
//
//    File: TRIGGER_init.cc
// Created: Wed Oct 24 06:29:48 EDT 2012
// Creator: davidl (on Darwin eleanor.jlab.org 12.2.0 i386)
//

#include <JANA/JEventLoop.h>
using namespace jana;

#include "DMCTrigger_factory.h"

jerror_t TRIGGER_init(JEventLoop *loop) {

	loop->AddFactory(new DMCTrigger_factory());

	return NOERROR;
}



