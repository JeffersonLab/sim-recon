// $Id$
//
//    File: DTrackFitterALT1.h
// Created: Tue Sep  2 11:18:22 EDT 2008
// Creator: davidl
//

#ifndef _DTrackFitterALT1_
#define _DTrackFitterALT1_

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
#include "DTrackFitter.h"

class DTrackCandidate;
class DTrack;
class DCDCTrackHit;
class DFDCPseudo;
class DMCThrown;

class DTrackFitterALT1:public DTrackFitter{
	public:
		DTrackFitterALT1(JEventLoop *loop);
		~DTrackFitterALT1();
	
		typedef DReferenceTrajectory::swim_step_t swim_step_t;
		
		// Virtual methods from TrackFitter base class
		string Name(void) const {return string("ALT1");}
		fit_status_t FitTrack(void);

	private:
		enum state_types{
			state_px,		///< x-momentum in RT coordinate system in GeV/c
			state_py,		///< y-momentum in RT coordinate system in GeV/c
			state_pz,		///< z-momentum in RT coordinate system in GeV/c
			state_x,			///< x-coordinate in RT coordinate system in cm
			state_v,			///< position-coordinate in RT coordinate system in cm perpendicular to x both and momentum direction
		};
		
		class hitInfo{
			public:
				const DCoordinateSystem* wire;	// Wire definitions
				double dist;							// Effective wire shifts due to drift time
				double err;								// Errors on drift time (or wire position) measurement
				double u_dist;							// Distances along the wire (for FDC cathodes)
				double u_err;							// Errors on distance along the wire (for FDC cathodes)
		};
		
		typedef vector<hitInfo> hitsInfo;

		double GetDistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double &s);
		double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL);
		double ChiSq(DReferenceTrajectory *rt, hitsInfo &hinfo, vector<double> &chisqv, double *chisq_ptr=NULL, int *dof_ptr=NULL);
		double ChiSq(DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt, hitsInfo &hinfo, vector<double> &chisqv, double *chisq_ptr=NULL, int *dof_ptr=NULL);
		void GetWiresShiftsErrs(fit_type_t fit_type, DReferenceTrajectory *rt, hitsInfo &hinfo);
		fit_status_t LeastSquaresB(hitsInfo &hinfo, DReferenceTrajectory *rt);
		void FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom);

		std::vector<double> chisqv;
		std::vector<double> sigmav;
		double Ngood_chisq_hits;
		DCoordinateSystem *target;
		DMatrix last_covariance;
		
		int eventnumber;
		const JGeometry *dgeom;
		double MAX_HIT_DIST;
		double CDC_Z_MIN;
		double CDC_Z_MAX;
		vector<DReferenceTrajectory*>rtv;
		DReferenceTrajectory *rt, *tmprt;
		bool hit_based;
		bool DEBUG_HISTS;
		int  DEBUG_LEVEL;
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
		bool TARGET_CONSTRAINT;
		
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
		TH1F *cdc_can_resi, *fdc_can_resi, *fdc_can_resi_cath;
		TH2F *chisq_vs_p_vs_theta;
};

#endif // _DTrackFitterALT1_

