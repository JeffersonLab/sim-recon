// $Id$
//
//    File: DNeutralParticleHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralParticleHypothesis_factory_
#define _DNeutralParticleHypothesis_factory_

#include <limits>

#include <JANA/JFactory.h>
#include <DANA/DApplication.h>
#include <HDGEOMETRY/DGeometry.h>
#include <PID/DNeutralParticleHypothesis.h>
#include <PID/DNeutralShower.h>
#include <PID/DEventRFBunch.h>
#include <PID/DParticleID.h>
#include <PID/DVertex.h>
#include <DVector3.h>
#include <DMatrixDSym.h>
#include <DMatrix.h>
#include <TMath.h>

class DNeutralParticleHypothesis_factory:public jana::JFactory<DNeutralParticleHypothesis>
{
	public:
		DNeutralParticleHypothesis_factory(){};
		~DNeutralParticleHypothesis_factory(){};

		DNeutralParticleHypothesis* Create_DNeutralParticleHypothesis(const DNeutralShower* locNeutralShower, Particle_t locPID, const DEventRFBunch* locEventRFBunch, const DVertex* locVertex) const;

		void Calc_ParticleCovariance_Photon(const DNeutralShower* locNeutralShower, const DVector3& locMomentum, const DVector3& locPathVector, DMatrixDSym& locParticleCovariance) const;
		void Calc_ParticleCovariance_Massive(const DNeutralShower* locNeutralShower, double locMass, double locDeltaT, double locStartTimeVariance, const DVector3& locMomentum, const DVector3& locPathVector, DMatrixDSym& locParticleCovariance) const;

	private:
		double dTargetLength;
		double dTargetRadius;
		DVector3 dTargetCenter;
		const DParticleID* dParticleID;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DNeutralParticleHypothesis_factory_

