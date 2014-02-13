#ifndef _DFactoryGenerator_monitoring_hists_
#define _DFactoryGenerator_monitoring_hists_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_monitoring_hists.h"

class DFactoryGenerator_monitoring_hists: public jana::JFactoryGenerator
{
	public:
		DFactoryGenerator_monitoring_hists(){}
		virtual ~DFactoryGenerator_monitoring_hists(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_monitoring_hists";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop)
		{
			loop->AddFactory(new DReaction_factory_monitoring_hists());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_monitoring_hists_

