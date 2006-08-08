// $Id$
//
//    File: DTrack_factory.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_factory_
#define _DTrack_factory_

#include <vector>

#include <TVector3.h>
#include <TMatrixD.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"

class DTrackCandidate;
class DTrack;
class DTrackHit;

class DTrack_factory:public JFactory<DTrack>{
	public:
		DTrack_factory(){};
		~DTrack_factory(){};
		const string toString(void);
	
	private:
		enum state_types{
			state_p,
			state_theta,
			state_phi,
			state_x,
			state_y
		};

		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t fini(void);

		typedef DReferenceTrajectory::swim_step_t swim_step_t;

		DTrack* FitTrack(const DTrackCandidate *trackcandidate);
		void GetTrackHits(DReferenceTrajectory *rt);
		TVector3 GetDistToRT(TVector3 &hit, swim_step_t *s2);
		void KalmanFilter(TMatrixD &state, TMatrixD &P, DReferenceTrajectory *rt);
		void KalmanStep(TMatrixD &x, TMatrixD &P_prev, TMatrixD &A, TMatrixD &H, TMatrixD &Q, TMatrixD &R, TMatrixD &W, TMatrixD &V);

		std::vector<const DTrackHit* > trackhits;
		std::vector<const DTrackHit* > hits_on_track;
		std::vector<TVector3> hit_dists;

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		string TRACKHIT_SOURCE;
		double MAX_HIT_DIST;
		DReferenceTrajectory::swim_step_t *swim_steps;
		int max_swim_steps;
};

#endif // _DTrack_factory_

