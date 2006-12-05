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
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"
#include "DCoordinateSystem.h"

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;
class DFDCPseudo;

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
		}cdc_hit_on_track_t;

		typedef struct{
			const DFDCPseudo* fdchit;
			const swim_step_t *swim_step;
			double dist_to_rt2;
			double dist;
			double s;
		}fdc_hit_on_track_t;

	private:
		enum state_types{
			state_px,		// x-momentum in RT coordinate system in GeV/c
			state_py,		// y-momentum in RT coordinate system in GeV/c
			state_pz,		// z-momentum in RT coordinate system in GeV/c
			state_x,			// x-coordinate in RT coordinate system in cm
			state_y,			// y-coordinate in RT coordinate system in cm
			state_z			// z-coordinate in RT coordinate system in cm
		};

		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t fini(void);

		DTrack* FitTrack(const DTrackCandidate *tc);
		void GetCDCTrackHits(DReferenceTrajectory *rt, double max_hit_dist=0.0);
		void GetFDCTrackHits(DReferenceTrajectory *rt, double max_hit_dist=0.0);
		double GetDistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double &s);
		//void KalmanFilter(TMatrixD &state, TMatrixD &P, DReferenceTrajectory *rt, TVector3 &vertex_pos, TVector3 &vertex_mom);
		//void KalmanStep(TMatrixD &x, TMatrixD &P, TMatrixD &z_minus_h, TMatrixD &A, TMatrixD &H, TMatrixD &Q, TMatrixD &R, TMatrixD &W, TMatrixD &V);
		//double ProjectStateBackwards(double q, TMatrixD &state, const swim_step_t *step, const swim_step_t *stepRT, const DCDCWire *wire);
		double ChiSq(double q, TMatrixD &state, swim_step_t *start_step, DReferenceTrajectory *rt=NULL);
		double ChiSq(double q, const TVector3 &pos, const TVector3 &mom, DReferenceTrajectory *rt=NULL);
		double LeastSquares(TVector3 &pos, TVector3 &mom, DReferenceTrajectory *rt, TVector3 &vertex_pos, TVector3 &vertex_mom);
		void FillDebugHists(DReferenceTrajectory *rt, TVector3 &vertex_pos, TVector3 &vertex_mom);

		std::vector<const DCDCTrackHit* > cdctrackhits;
		std::vector<const DFDCPseudo* > fdctrackhits;
		std::vector<cdc_hit_on_track_t > cdchits_on_track;
		std::vector<fdc_hit_on_track_t > fdchits_on_track;
		
		std::vector<double> chisqv;
		std::vector<double> sigmav;

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		string TRACKHIT_SOURCE;
		double MAX_HIT_DIST;
		double CDC_Z_MIN;
		double CDC_Z_MAX;
		swim_step_t *swim_steps;
		int max_swim_steps;
		swim_step_t *swim_steps_ls;
		int max_swim_steps_ls;
		bool hit_based;
		bool DEBUG_HISTS;
		
		TH3F *cdcdoca_vs_dist_vs_ring;
		TH2F *cdcdoca_vs_dist;
		TH2F *fdcdoca_vs_dist;
		TH1F *dist_stereo, *dist_axial;
		TH1F *doca_stereo, *doca_axial;
		TH1F *fit_chisq;
		TH1F *residuals;
};

#endif // _DTrack_factory_

