#ifndef _DFactoryGenerator_DReaction_
#define _DFactoryGenerator_DReaction_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory.h"

class DFactoryGenerator_DReaction: public jana::JFactoryGenerator
{
	public:
		DFactoryGenerator_DReaction(){}
		virtual ~DFactoryGenerator_DReaction(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_DReaction";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop)
		{
			loop->AddFactory(new DReaction_factory());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_DReaction_

