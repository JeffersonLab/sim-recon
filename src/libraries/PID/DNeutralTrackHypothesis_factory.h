// $Id$
//
//    File: DNeutralTrackHypothesis_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralTrackHypothesis_factory_
#define _DNeutralTrackHypothesis_factory_

#include <JANA/JFactory.h>
#include <PID/DNeutralTrackHypothesis.h>
#include <PID/DParticleID.h>
#include <PID/DNeutralShowerCandidate.h>
#include <PID/DChargedTrackHypothesis.h>
#include <PID/DVertex.h>
#include <PID/DChargedTrack.h>
#include <PID/DKinematicData.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <DVector3.h>
#include <DMatrixDSym.h>

class DNeutralTrackHypothesis_factory:public jana::JFactory<DNeutralTrackHypothesis>{
	public:
		DNeutralTrackHypothesis_factory(){};
		~DNeutralTrackHypothesis_factory(){};

		void Build_ErrorMatrix(const DVector3 &locPathVector, double locEnergy, const DMatrixDSym& locVariances, DMatrixDSym& locErrorMatrix);
		void Calc_Variances(const DNeutralShowerCandidate *locNeutralShowerCandidate, double locParticleEnergyUncertainty, DMatrixDSym &locVariances);

	private:
		DParticleID *dPIDAlgorithm;
		double dTargetLength;
		double dTargetRadius;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DNeutralTrackHypothesis_factory_

