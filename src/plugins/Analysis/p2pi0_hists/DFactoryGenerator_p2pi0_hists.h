// $Id$
//
//    File: DFactoryGenerator_p2pi0_hists.h
// Created: Tue May 19 07:56:15 EDT 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_p2pi0_hists_
#define _DFactoryGenerator_p2pi0_hists_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_p2pi0_hists.h"

class DFactoryGenerator_p2pi0_hists : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_p2pi0_hists";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_p2pi0_hists());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_p2pi0_hists_

