// $Id$
//
//    File: DFactory_DHDDMBCALTruth.h
// Created: Fri Nov 18 10:37:37 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DFactory_DHDDMBCALTruth_
#define _DFactory_DHDDMBCALTruth_

#include "DFactory.h"
#include "DHDDMBCALTruth.h"

class DFactory_DHDDMBCALTruth:public DFactory<DHDDMBCALTruth>{
	public:
		DFactory_DHDDMBCALTruth(){};
		~DFactory_DHDDMBCALTruth(){};

                derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DHDDMBCALTruth_

