// $Id$
//
//    File: DFactoryGenerator_ee_convert.h
// Created: Wed Jun 14 06:17:54 EDT 2017
// Creator: jrsteven (on Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#ifndef _DFactoryGenerator_ee_convert_
#define _DFactoryGenerator_ee_convert_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_ee_convert.h"

class DFactoryGenerator_ee_convert : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_ee_convert";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_ee_convert());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_ee_convert_

