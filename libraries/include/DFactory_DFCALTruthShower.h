// $Id$
//
//    File: DFactory_DFCALTruthShower.h
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#ifndef _DFactory_DFCALTruthShower_
#define _DFactory_DFCALTruthShower_

#include "DFactory.h"
#include "DFCALTruthShower.h"

class DFactory_DFCALTruthShower:public DFactory<DFCALTruthShower>{
	public:
		DFactory_DFCALTruthShower(){};
		~DFactory_DFCALTruthShower(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);


	private:
		//derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DFCALTruthShower_

