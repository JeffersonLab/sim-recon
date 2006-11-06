// $Id$
//
//    File: DUPVTruthHit_factory.h
// Created: Mon Nov  6 09:58:35 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.8.0 powerpc)
//

#ifndef _DUPVTruthHit_factory_
#define _DUPVTruthHit_factory_

#include <JANA/JFactory.h>
#include "DUPVTruthHit.h"

class DUPVTruthHit_factory:public JFactory<DUPVTruthHit>{
	public:
		DUPVTruthHit_factory(){};
		~DUPVTruthHit_factory(){};
		const string toString(void);

		double GetETotal(void);

	private:
};

#endif // _DUPVTruthHit_factory_

