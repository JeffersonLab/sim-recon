// $Id$
//
//    File: DFactoryGenerator_ppi0g_hists.h
//

#ifndef _DFactoryGenerator_ppi0g_hists_
#define _DFactoryGenerator_ppi0g_hists_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_ppi0g_hists.h"

class DFactoryGenerator_ppi0g_hists : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_ppi0g_hists";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_ppi0g_hists());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_ppi0g_hists_

