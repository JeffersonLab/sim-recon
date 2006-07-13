// $Id$
//
//    File: DFCALTruthShower_factory.h
// Created: Wed Jan  4 14:43:05 EST 2006
// Creator: davidl (on Linux jlabl1.jlab.org 2.4.21-37.ELsmp i686)
//

#ifndef _DFCALTruthShower_factory_
#define _DFCALTruthShower_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"
#include "DFCALTruthShower.h"

class DFCALTruthShower_factory:public JFactory<DFCALTruthShower>{
	public:
		DFCALTruthShower_factory(){};
		~DFCALTruthShower_factory(){};
		jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);


	private:
		//jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DFCALTruthShower_factory_

