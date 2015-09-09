// $Id$
//
//    File: DMCThrownMatching_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DMCThrownMatching_factory_
#define _DMCThrownMatching_factory_

#include <map>

#include "TDecompLU.h"

#include <JANA/JFactory.h>
#include "TRACKING/DMCThrown.h"
#include "PID/DChargedTrackHypothesis.h"
#include "PID/DChargedTrack.h"
#include "PID/DNeutralParticleHypothesis.h"
#include "PID/DNeutralParticle.h"
#include "PID/DNeutralShower.h"
#include "PID/DBeamPhoton.h"
#include "ANALYSIS/DMCThrownMatching.h"
#include "DANA/DApplication.h"

#include "TRACKING/DTrackTimeBased.h"
#include "TAGGER/DTAGMHit.h"
#include "TAGGER/DTAGHHit.h"
#include "TOF/DTOFPoint.h"
#include "TOF/DTOFTruth.h"
#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALTruthShower.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALTruthShower.h"

using namespace jana;
using namespace std;

class DMCThrownMatching_factory : public jana::JFactory<DMCThrownMatching>
{
	public:
		bool Calc_InverseMatrix(const DMatrixDSym& locInputCovarianceMatrix, DMatrixDSym& locInverse3x3Matrix) const;
		double Calc_MatchFOM(const DVector3& locMomentum_Thrown, const DVector3& locMomentum_Detected, DMatrixDSym locInverse3x3Matrix) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void Find_GenReconMatches_BeamPhotons(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const;

		void Find_GenReconMatches_ChargedTrack(const vector<const DChargedTrack*>& locChargedTracks, DMCThrownMatching* locMCThrownMatching) const;
		void Find_GenReconMatches_ChargedHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DChargedTrackHypothesis*>& locInputChargedTrackHypothesisVector, DMCThrownMatching* locMCThrownMatching) const;

		void Find_GenReconMatches_NeutralParticle(const vector<const DNeutralParticle*>& locNeutralParticles, DMCThrownMatching* locMCThrownMatching) const;
		void Find_GenReconMatches_NeutralHypo(const vector<const DMCThrown*>& locInputMCThrownVector, const vector<const DNeutralParticleHypothesis*>& locInputNeutralParticleHypothesisVector, DMCThrownMatching* locMCThrownMatching) const;

		void Find_GenReconMatches_TOFPoints(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const;
		void Find_GenReconMatches_BCALShowers(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const;
		void Find_GenReconMatches_FCALShowers(JEventLoop* locEventLoop, DMCThrownMatching* locMCThrownMatching) const;

		double dMinTrackMatchHitFraction;
		double dMaximumTOFMatchDistance;
		double dMaximumFCALMatchDistance;
		double dMaximumBCALMatchAngleDegrees;
		unsigned int dDebugLevel;
		double dTargetCenter;
		double dMaxTotalParticleErrorForMatch;
};

#endif // _DMCThrownMatching_factory_

