// $Id$
//
//    File: DNeutralParticle_factory_PreSelect.h
// Created: Mon Dec  7 14:29:24 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticle_factory_PreSelect_
#define _DNeutralParticle_factory_PreSelect_

#include <JANA/JFactory.h>
#include <PID/DNeutralParticle.h>
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DNeutralShower.h>

using namespace std;
using namespace jana;

class DNeutralParticle_factory_PreSelect : public jana::JFactory<DNeutralParticle>
{
	public:
		DNeutralParticle_factory_PreSelect(){};
		~DNeutralParticle_factory_PreSelect(){};
		const char* Tag(void){return "PreSelect";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DNeutralParticle_factory_PreSelect_

