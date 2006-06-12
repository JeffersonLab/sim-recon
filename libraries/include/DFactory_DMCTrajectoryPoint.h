// $Id$
//
//    File: DFactory_DMCTrajectoryPoint.h
// Created: Mon Jun 12 09:29:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.6.0 powerpc)
//

#ifndef _DFactory_DMCTrajectoryPoint_
#define _DFactory_DMCTrajectoryPoint_

#include "DFactory.h"
#include "DMCTrajectoryPoint.h"

class DFactory_DMCTrajectoryPoint:public DFactory<DMCTrajectoryPoint>{
	public:
		DFactory_DMCTrajectoryPoint(){};
		~DFactory_DMCTrajectoryPoint(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
};

#endif // _DFactory_DMCTrajectoryPoint_

