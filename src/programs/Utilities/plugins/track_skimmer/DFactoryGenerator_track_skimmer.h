// $Id$
//
//    File: DFactoryGenerator_track_skimmer.h
// Created: Tue Jan 13 11:08:16 EST 2015
// Creator: Paul (on Darwin Pauls-MacBook-Pro.local 14.0.0 i386)
//

#ifndef _DFactoryGenerator_track_skimmer_
#define _DFactoryGenerator_track_skimmer_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_track_skimmer.h"

class DFactoryGenerator_track_skimmer : public jana::JFactoryGenerator
{
	public:
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_track_skimmer";}
		
		jerror_t GenerateFactories(jana::JEventLoop* locEventLoop)
		{
			locEventLoop->AddFactory(new DReaction_factory_track_skimmer());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_track_skimmer_

