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
#include "CDC/DCDCWire.h"

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;

class DTrack_factory:public JFactory<DTrack>{
	public:
		DTrack_factory(){};
		~DTrack_factory(){};
		const string toString(void);
	
		typedef DReferenceTrajectory::swim_step_t swim_step_t;

		typedef struct{
			const DCDCTrackHit* cdchit;
			const swim_step_t *swim_step;
			double dist_to_rt2;
			double dist;
			double s;
		}hit_on_track_t;

	private:
		enum state_types{
			state_p,			// Total momentum in GeV/c
			state_theta,	// theta-angle of momentum in lab coordinate system
			state_phi,		// phi-angle of momentum in lab coordinate system
			state_x,			// x-coordinate in RT coordinate system
			state_y			// y-coordinate in RT coordinate system
		};

		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t fini(void);

		DTrack* FitTrack(const DTrackCandidate *trackcandidate);
		void GetCDCTrackHits(DReferenceTrajectory *rt);
		double GetDistToRT(const DCDCWire *wire, const swim_step_t *step, double &s);
		swim_step_t* KalmanFilter(TMatrixD &state, TMatrixD &P, DReferenceTrajectory *rt);
		void KalmanStep(TMatrixD &x, TMatrixD &P, TMatrixD &z_minus_h, TMatrixD &A, TMatrixD &H, TMatrixD &Q, TMatrixD &R, TMatrixD &W, TMatrixD &V);
		double ProjectState(double q, TMatrixD &state, const swim_step_t *step, const swim_step_t *stepRT, const DCDCWire *wire);

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

