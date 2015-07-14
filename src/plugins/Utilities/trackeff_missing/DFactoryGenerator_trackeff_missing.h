// $Id$
//
//    File: DFactoryGenerator_trackeff_missing.h
// Created: Wed Feb 25 08:58:19 EST 2015
// Creator: pmatt (on Linux pmattdesktop.jlab.org 2.6.32-504.8.1.el6.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_trackeff_missing_
#define _DFactoryGenerator_trackeff_missing_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_trackeff_missing.h"

class DFactoryGenerator_trackeff_missing : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_trackeff_missing";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_trackeff_missing());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_trackeff_missing_

