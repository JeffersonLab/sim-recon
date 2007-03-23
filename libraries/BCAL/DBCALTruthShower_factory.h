// $Id$
//
//    File: DBCALTruthShower_factory.h
// Created: Fri Nov 18 10:37:37 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALTruthShower_factory_
#define _DBCALTruthShower_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"
#include "DBCALTruthShower.h"

class DBCALTruthShower_factory:public JFactory<DBCALTruthShower>{
	public:
		DBCALTruthShower_factory(){};
		~DBCALTruthShower_factory(){};

		const string toString(void);
};

#endif // _DBCALTruthShower_factory_

