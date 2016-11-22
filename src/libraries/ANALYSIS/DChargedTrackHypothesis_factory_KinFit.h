// $Id$
//
//    File: DChargedTrackHypothesis_factory_KinFit.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DChargedTrackHypothesis_factory_KinFit_
#define _DChargedTrackHypothesis_factory_KinFit_

#include <JANA/JFactory.h>
#include <PID/DChargedTrackHypothesis.h>
#include <PID/DChargedTrack.h>
#include <KINFITTER/DKinFitParticle.h>
#include <ANALYSIS/DParticleCombo.h>
#include <PID/DParticleID.h>

using namespace jana;
using namespace std;

class DChargedTrackHypothesis_factory_KinFit : public jana::JFactory<DChargedTrackHypothesis>
{
	public:
		DChargedTrackHypothesis_factory_KinFit(){};
		~DChargedTrackHypothesis_factory_KinFit(){};
		const char* Tag(void){return "KinFit";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DChargedTrackHypothesis* Build_ChargedTrackHypothesis(const DChargedTrackHypothesis* locChargedTrackHypothesis, DKinFitParticle* locKinFitParticle, const DChargedTrack* locChargedTrack, const DParticleCombo* locParticleCombo);

		const DParticleID* dPIDAlgorithm;
};

#endif // _DChargedTrackHypothesis_factory_KinFit_

