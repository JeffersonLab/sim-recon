// $Id$
//
//    File: DTrack_factory.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_factory_
#define _DTrack_factory_

#include <vector>

#include <DVector3.h>
#include <DMatrix.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JGeometry.h>
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "DTrack.h"
#include "DReferenceTrajectory.h"
#include "DCoordinateSystem.h"

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;
class DFDCPseudo;
class DMCThrown;

class DTrack_factory:public JFactory<DTrack>{
	public:
		DTrack_factory();
		~DTrack_factory();
		const string toString(void);
	
		typedef DReferenceTrajectory::swim_step_t swim_step_t;
#if 0
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
#endif
	private:
		enum state_types{
			state_px,		///< x-momentum in RT coordinate system in GeV/c
			state_py,		///< y-momentum in RT coordinate system in GeV/c
			state_pz,		///< z-momentum in RT coordinate system in GeV/c
			state_x,			///< x-coordinate in RT coordinate system in cm
			state_v,			///< position-coordinate in RT coordinate system in cm perpendicular to x both and momentum direction
		};
		
		enum fit_status_t{
			FIT_OK,
			FIT_FAILED
		};

		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t fini(void);

		void FindHitCandidateProbabilities(void);
		DTrack* FitTrack(DReferenceTrajectory* rt, int candidateid);
		void GetCDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob);
		void GetFDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob);
		double GetDistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double &s);
		//void KalmanFilter(DMatrix &state, DMatrix &P, DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom);
		//void KalmanStep(DMatrix &x, DMatrix &P, DMatrix &z_minus_h, DMatrix &A, DMatrix &H, DMatrix &Q, DMatrix &R, DMatrix &W, DMatrix &V);
		//double ProjectStateBackwards(double q, DMatrix &state, const swim_step_t *step, const swim_step_t *stepRT, const DCDCWire *wire);
		double ChiSq(double q, DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt=NULL);
		double ChiSq(double q, const DVector3 &pos, const DVector3 &mom, DReferenceTrajectory *rt=NULL);
		double ChiSq(double q, DReferenceTrajectory *rt);
		fit_status_t LeastSquares(DVector3 &pos, DVector3 &mom, DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom, double &chisq);
		void FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom);

		std::vector<const DTrackCandidate*>	trackcandidates;
		std::vector<const DCDCTrackHit* >	cdctrackhits;
		std::vector<const DFDCPseudo* >		fdctrackhits;
		std::vector<const DCDCTrackHit* >	cdchits_on_track;
		std::vector<const DFDCPseudo* >		fdchits_on_track;
		std::vector<vector<double> > cdcprobs;
		std::vector<vector<double> > fdcprobs;
		
		std::vector<double> chisqv;
		std::vector<double> sigmav;
		double Ngood_chisq_hits;
		DCoordinateSystem *target;
		DMatrix last_covariance;

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		double MAX_HIT_DIST;
		double CDC_Z_MIN;
		double CDC_Z_MAX;
		vector<DReferenceTrajectory*>rtv;
		DReferenceTrajectory *tmprt;
		//swim_step_t *swim_steps;
		//int max_swim_steps;
		//swim_step_t *swim_steps_ls;
		//int max_swim_steps_ls;
		bool hit_based;
		bool DEBUG_HISTS;
		bool USE_CDC;
		bool USE_FDC_ANODE;
		bool USE_FDC_CATHODE;
		double MAX_CHISQ_DIFF;
		int MAX_FIT_ITERATIONS;
		double SIGMA_CDC;
		double SIGMA_FDC_ANODE;
		double SIGMA_FDC_CATHODE;
		double CHISQ_MAX_RESI_SIGMAS;
		double CHISQ_GOOD_LIMIT;
		double LEAST_SQUARES_DP;
		double LEAST_SQUARES_DX;
		unsigned int LEAST_SQUARES_MIN_HITS;
		double LEAST_SQUARES_MAX_E2NORM;
		double DEFAULT_STEP_SIZE;
		double MIN_CDC_HIT_PROB;
		double MAX_CDC_DOUBLE_HIT_PROB;
		double MIN_FDC_HIT_PROB;
		double MAX_FDC_DOUBLE_HIT_PROB;
		double TOF_MASS;
		
		TH3F *cdcdoca_vs_dist_vs_ring;
		TH2F *cdcdoca_vs_dist;
		TH2F *fdcdoca_vs_dist, *fdcu_vs_s;
		TH1F *dist_stereo, *dist_axial;
		TH1F *doca_stereo, *doca_axial;
		TH2F *chisq_final_vs_initial;
		TH2F *nhits_final_vs_initial;
		TH1F *Npasses;
		TH1F *ptotal;
		TH2F *residuals_cdc, *residuals_fdc_anode, *residuals_fdc_cathode;
		TH2F *initial_chisq_vs_Npasses, *chisq_vs_pass, *dchisq_vs_pass;
		TH3F *residuals_cdc_vs_s, *residuals_fdc_anode_vs_s, *residuals_fdc_cathode_vs_s;
		TH1F *cdc_single_hit_prob, *cdc_double_hit_prob;
		TH1F *fdc_single_hit_prob, *fdc_double_hit_prob;
};

#endif // _DTrack_factory_

