//
//    File: DTwoGammaFit_factory.h
// Created: Thu Sep 13 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_factory_ETA_
#define _DTwoGammaFit_factory_ETA_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTwoGammaFit_factory.h"

class DTwoGammaFit_factory_ETA:public DTwoGammaFit_factory {
	public:
		DTwoGammaFit_factory_ETA();
		~DTwoGammaFit_factory_ETA(){};
		const char* Tag(void){return "ETA";}
	
	private:


};

#endif // _DTwoGammaFit_factory_

