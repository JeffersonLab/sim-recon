// $Id$
//
//    File: DFactoryGenerator_p2pi.h
// Created: Thu May  7 16:22:04 EDT 2015
// Creator: mstaib (on Linux gluon109.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_p2pi_
#define _DFactoryGenerator_p2pi_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_p2pi.h"

class DFactoryGenerator_p2pi : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_p2pi";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_p2pi());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_p2pi_

