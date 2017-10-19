// $Id$
//
//    File: DBeamPhoton_factory_TAGGEDMCGEN.h
// Created: Mon Aug  5 14:29:24 EST 2014
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DBeamPhoton_factory_TAGGEDMCGEN_
#define _DBeamPhoton_factory_TAGGEDMCGEN_

#include <JANA/JFactory.h>
#include <PID/DBeamPhoton.h>

class DBeamPhoton_factory_TAGGEDMCGEN:public jana::JFactory<DBeamPhoton>{
	public:
		const char* Tag(void){return "TAGGEDMCGEN";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
};

#endif // _DBeamPhoton_factory_TAGGEDMCGEN_

