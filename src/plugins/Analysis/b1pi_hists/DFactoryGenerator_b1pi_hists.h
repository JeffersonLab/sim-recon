#ifndef _DFactoryGenerator_b1pi_hists_
#define _DFactoryGenerator_b1pi_hists_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>

#include "DReaction_factory_b1pi_hists.h"
#include "DEventWriterROOT_factory_test.h"

class DFactoryGenerator_b1pi_hists: public jana::JFactoryGenerator
{
	public:
		DFactoryGenerator_b1pi_hists(){}
		virtual ~DFactoryGenerator_b1pi_hists(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGenerator_b1pi_hists";}
		
		jerror_t GenerateFactories(jana::JEventLoop *loop)
		{
			loop->AddFactory(new DReaction_factory_b1pi_hists());
	        loop->AddFactory(new DEventWriterROOT_factory_test());
			return NOERROR;
		}
};

#endif // _DFactoryGenerator_b1pi_hists_

