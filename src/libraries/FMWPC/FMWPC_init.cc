// $Id$
//
//    File: FMWPC_init.cc
// Created: Tue Jun 16 07:04:58 EDT 2015
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <JANA/JEventLoop.h>
#include <JANA/JFactory.h>
using namespace jana;

#include "DFMWPCHit.h"
#include "DFMWPCTruthHit.h"

jerror_t FMWPC_init(JEventLoop *loop) {

	/// Create and register FMWPC data factories
	loop->AddFactory(new JFactory<DFMWPCHit>());
	loop->AddFactory(new JFactory<DFMWPCTruthHit>());

	return NOERROR;
}

