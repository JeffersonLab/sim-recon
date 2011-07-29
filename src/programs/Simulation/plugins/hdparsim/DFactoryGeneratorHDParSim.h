// $Id$
//
//    File: DFactoryGeneratorHDParSim.h
// Created: Tue Feb  3 11:25:33 EST 2009
// Creator: davidl (on Darwin Harriet.local 9.6.0 i386)
//

#ifndef _DFactoryGeneratorHDParSim_
#define _DFactoryGeneratorHDParSim_

#include <JANA/jerror.h>
#include <JANA/JFactoryGenerator.h>
using namespace jana;

#include "DTrackTimeBased_factory_HDParSim.h"
#include "DPhoton_factory_HDParSim.h"

class DFactoryGeneratorHDParSim: public JFactoryGenerator{
	public:
		DFactoryGeneratorHDParSim(){pthread_mutex_init(&root_mutex, NULL);}
		virtual ~DFactoryGeneratorHDParSim(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DFactoryGeneratorHDParSim";}
		
		jerror_t GenerateFactories(JEventLoop *loop){
			pthread_mutex_lock(&root_mutex);
			
			loop->AddFactory(new DTrackTimeBased_factory_HDParSim());
			loop->AddFactory(new DPhoton_factory_HDParSim());

			pthread_mutex_unlock(&root_mutex);

			return NOERROR;
		}

	protected:
	
	
	private:
		pthread_mutex_t root_mutex;

};

#endif // _DFactoryGeneratorHDParSim_

