// $Id$
//
//    File: DNeutralParticleHypothesis_factory_KinFit.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticleHypothesis_factory_KinFit_
#define _DNeutralParticleHypothesis_factory_KinFit_

#include <deque>

#include "JANA/JFactory.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralShower.h"

#include "ANALYSIS/DKinFitParticle.h"
#include "ANALYSIS/DParticleCombo.h"

using namespace jana;
using namespace std;

class DNeutralParticleHypothesis_factory_KinFit : public jana::JFactory<DNeutralParticleHypothesis>
{
	public:
		DNeutralParticleHypothesis_factory_KinFit(){};
		~DNeutralParticleHypothesis_factory_KinFit(){};
		const char* Tag(void){return "KinFit";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DNeutralParticleHypothesis* Build_NeutralParticleHypothesis(const DNeutralParticleHypothesis* locNeutralParticleHypothesis, const DKinFitParticle* locKinFitParticle, const DNeutralShower* locNeutralShower, const DParticleCombo* locParticleCombo);
};

#endif // _DNeutralParticleHypothesis_factory_KinFit_

