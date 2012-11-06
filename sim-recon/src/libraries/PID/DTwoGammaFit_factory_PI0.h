//
//    File: DTwoGammaFit_factory.h
// Created: Tue Apr 17 11:57:50 EST 2007
// Creator: M. Kornicer (on Linux stan)
//

#ifndef _DTwoGammaFit_factory_PI0_
#define _DTwoGammaFit_factory_PI0_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTwoGammaFit_factory.h"

class DTwoGammaFit_factory_PI0:public DTwoGammaFit_factory {
	public:
		DTwoGammaFit_factory_PI0();
		~DTwoGammaFit_factory_PI0(){};
      const char* Tag(void){return "PI0";}
	
	private:


};

#endif // _DTwoGammaFit_factory_

