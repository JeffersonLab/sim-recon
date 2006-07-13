// $Id$
//
//    File: DHDDMBCALTruth_factory.h
// Created: Fri Nov 18 10:37:37 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DHDDMBCALTruth_factory_
#define _DHDDMBCALTruth_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"
#include "DHDDMBCALTruth.h"

class DHDDMBCALTruth_factory:public JFactory<DHDDMBCALTruth>{
	public:
		DHDDMBCALTruth_factory(){};
		~DHDDMBCALTruth_factory(){};

                jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);

		const string toString(void);


	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DHDDMBCALTruth_factory_

