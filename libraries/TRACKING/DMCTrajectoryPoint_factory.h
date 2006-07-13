// $Id$
//
//    File: DMCTrajectoryPoint_factory.h
// Created: Mon Jun 12 09:29:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.6.0 powerpc)
//

#ifndef _DMCTrajectoryPoint_factory_
#define _DMCTrajectoryPoint_factory_

#include "JANA/JFactory.h"
#include "DMCTrajectoryPoint.h"

class DMCTrajectoryPoint_factory:public JFactory<DMCTrajectoryPoint>{
	public:
		DMCTrajectoryPoint_factory(){};
		~DMCTrajectoryPoint_factory(){};
		const string toString(void);


	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
};

#endif // _DMCTrajectoryPoint_factory_

