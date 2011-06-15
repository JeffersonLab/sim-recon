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
#include "DTrackWireBased.h"
#include "DReferenceTrajectory.h"
#include "DCoordinateSystem.h"
#include "DTrackFitter.h"

class DTrackCandidate;
class DTrackWireBased;
class DCDCTrackHit;
class DFDCPseudo;
class DMCThrown;

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="ND.png" width="100">
///	<IMG src="DEP.png" width="100">
///	</A>
/// \endhtmlonly

/// A global least-squares track fitter. This has been superceded by the DTrackFitterKalmanSIMD class.

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
		
		enum resi_types{
			resi_type_cdc_anode,
			resi_type_fdc_anode,
			resi_type_fdc_cathode,
			resi_type_other
		};
		
		class hitInfo{
			public:
				const DCoordinateSystem* wire;	///< Wire definitions
				double dist;							///< Effective wire shifts due to drift time
				double err;								///< Errors on drift time (or wire position) measurement
				double u_dist;							///< Distances along the wire (for FDC cathodes)
				double u_lorentz;						///< Lorentz correction  to u_dist ( u = u_dist + u_lorentz )
				double u_err;							///< Errors on distance along the wire (for FDC cathodes)
				bool   good;							///< Set by GetResiInfo if dist is used
				bool   good_u;							///< Set by GetResiInfo if u_dist is used
		};
		
		typedef vector<hitInfo> hitsInfo;
		
		class resiInfo{
			public:
				hitInfo *hit;
				int layer;
				int resi_type;
				double resi;
				double err;
				const swim_step_t *step;
		};

		double ChiSq(fit_type_t fit_type, DReferenceTrajectory *rt, double *chisq_ptr=NULL, int *dof_ptr=NULL, vector<pull_t> *pulls_ptr=NULL);
		double ChiSq(vector<resiInfo> &residuals, double *chisq_ptr=NULL, int *dof_ptr=NULL);
		void GetHits(fit_type_t fit_type, DReferenceTrajectory *rt, hitsInfo &hinfo);
		vector<bool> GetResiInfo(DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt, hitsInfo &hinfo, vector<resiInfo> &residuals);
		vector<bool> GetResiInfo(DReferenceTrajectory *rt, hitsInfo &hinfo, vector<resiInfo> &residuals);
		fit_status_t LeastSquaresB(hitsInfo &hinfo, DReferenceTrajectory *rt);
		void FilterGood(DMatrix &my_resiv, vector<bool> &my_good, vector<bool> &good_all);
		void PrintChisqElements(DMatrix &resiv, DMatrix &cov_meas, DMatrix &cov_muls, DMatrix &weights);
		void ForceLRTruth(JEventLoop *loop, DReferenceTrajectory *rt, hitsInfo &hinfo);
		void FillDebugHists(DReferenceTrajectory *rt, DVector3 &vertex_pos, DVector3 &vertex_mom);

		// The following are filled by the last call to one of the
		// ChiSq(...) methods
		DMatrix resiv;		///< residuals vector (Nmeasurements x 1)
		DMatrix cov_meas;	///< Measurement errors of hits (diagonal Nmeasurements x Nmeasurements)
		DMatrix cov_muls;	///< Covariance of hits due to multiple scattering (Nmeasurements x Nmeasurements)
		DMatrix weights;	///< Inverse of cov_meas + cov_muls

		// The following are filled by the last call to LeastSquaresB(...)
		DMatrix cov_parm;	///< Covariance of fit parameters (Nparms x Nparms (where Nparms=5))

		// Other data members
		DCoordinateSystem *target;
		
		int eventnumber;
		const JGeometry *dgeom;
		vector<DReferenceTrajectory*>rtv;
		DReferenceTrajectory *rt, *tmprt;
		bool hit_based;
		bool DEBUG_HISTS;
		int  DEBUG_LEVEL;
		double MAX_CHISQ_DIFF;
		int MAX_FIT_ITERATIONS;
		double SIGMA_CDC;
		bool CDC_USE_PARAMETERIZED_SIGMA;
		double SIGMA_FDC_ANODE;
		double SIGMA_FDC_CATHODE;
		double CHISQ_GOOD_LIMIT;
		double LEAST_SQUARES_DP;
		double LEAST_SQUARES_DX;
		unsigned int LEAST_SQUARES_MIN_HITS;
		double LEAST_SQUARES_MAX_E2NORM;
		double DEFAULT_MASS;
		bool TARGET_CONSTRAINT;
		bool LR_FORCE_TRUTH;
		bool USE_MULS_COVARIANCE;
		bool USE_FDC;
		bool USE_FDC_CATHODE;
		bool USE_CDC;
		string MATERIAL_MAP_MODEL;
		
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
		TH1F *lambda;
};

#endif // _DTrackFitterALT1_

