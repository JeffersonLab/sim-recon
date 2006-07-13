// $Id$
//
//    File: DMCThrown_factory.h
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DMCThrown_factory_
#define _DMCThrown_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "HDDM/hddm_s.h"
#include "DMCThrown.h"

class DMCThrown_factory:public JFactory<DMCThrown>{
	public:
		DMCThrown_factory(){};
		~DMCThrown_factory(){};
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DMCThrown_factory_

