// $Id$
//
//    File: DFactoryGenerator_B3pi_eff_misspip.h
// Created: Fri Jun 30 00:38:13 EDT 2017
// Creator: jmhardin (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_B3pi_eff_misspip_
#define _DFactoryGenerator_B3pi_eff_misspip_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_B3pi_eff_misspip.h"

class DFactoryGenerator_B3pi_eff_misspip : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_B3pi_eff_misspip";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_B3pi_eff_misspip());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_B3pi_eff_misspip_

