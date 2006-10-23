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
#include <TH2.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;

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
		void GetCDCTrackHits(DReferenceTrajectory *rt);
		double GetDistToRT(const DCDCTrackHit *cdchit, const swim_step_t *step);
		void KalmanFilter(TMatrixD &state, TMatrixD &P, DReferenceTrajectory *rt);
		void KalmanStep(TMatrixD &x, TMatrixD &P, TMatrixD &z_minus_h, TMatrixD &A, TMatrixD &H, TMatrixD &Q, TMatrixD &R, TMatrixD &W, TMatrixD &V);

		typedef struct{
			const DCDCTrackHit* cdchit;
			const swim_step_t *swim_step;
			float dist_to_rt2;
		}hit_on_track_t;

		std::vector<const DCDCTrackHit* > cdctrackhits;
		std::vector<hit_on_track_t > cdchits_on_track;

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		string TRACKHIT_SOURCE;
		double MAX_HIT_DIST;
		double CDC_Z_MIN;
		double CDC_Z_MAX;
		swim_step_t *swim_steps;
		int max_swim_steps;
		
		TH1F *cdcdocart,*cdcdocaswim, *cdcdocatdrift;
};

#endif // _DTrack_factory_

