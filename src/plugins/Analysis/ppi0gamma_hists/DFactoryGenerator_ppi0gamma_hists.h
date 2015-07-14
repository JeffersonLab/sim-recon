// $Id$
//
//    File: DFactoryGenerator_ppi0gamma_hists.h
// Created: Fri May 15 14:19:50 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_ppi0gamma_hists_
#define _DFactoryGenerator_ppi0gamma_hists_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_ppi0gamma_hists.h"

class DFactoryGenerator_ppi0gamma_hists : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_ppi0gamma_hists";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_ppi0gamma_hists());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_ppi0gamma_hists_

