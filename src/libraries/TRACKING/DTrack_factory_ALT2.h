// $Id$
//
//    File: DTrack_factory_ALT2.h
// Created: Fri. March. 4, 2008
// Creator: davidl
//

#ifndef _DTrack_factory_ALT2_
#define _DTrack_factory_ALT2_

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

class DTrack_factory_ALT2:public JFactory<DTrack>{
	public:
		DTrack_factory_ALT2();
		~DTrack_factory_ALT2();
		const string toString(void);
		const char* Tag(void){return "ALT2";}
	
		typedef DReferenceTrajectory::swim_step_t swim_step_t;

	private:
		enum state_types{
			state_curvature,		///< 1/momentum in RT coordinate system in c/GeV
			state_dip,		///< dip angle of momentum in RT coordinate system in radians
			state_phi,		///< phi angle of momentum in RT coordinate system in radians
			state_x,			///< x-coordinate in RT coordinate system in cm
			state_v,			///< position-coordinate in RT coordinate system in cm perpendicular to x both and momentum direction
		};
		
		enum fit_status_t{
			FIT_OK,
			FIT_NO_IMPROVEMENT,
			FIT_FAILED
		};

		enum fit_type_t{
			kHitBased,
			kTimeBased
		};

		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
		jerror_t fini(void);

		void FindHitCandidateProbabilities(void);
		DTrack* FitTrack(DReferenceTrajectory* rt, int candidateid);
		DTrack* FitTrackWithOppositeCharge(DReferenceTrajectory* rt, int candidateid, DTrack* &track);
		void GetCDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob);
		void GetFDCTrackHitProbabilities(DReferenceTrajectory *rt, vector<double> &prob);
		double GetDistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double &s);
		double ChiSq(double q, DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt=NULL);
		double ChiSq(double q, const DVector3 &pos, const DVector3 &mom, DReferenceTrajectory *rt=NULL);
		double ChiSq(double q, DReferenceTrajectory *rt);
		double ChiSq(DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv);
		double ChiSq(DMatrix &state, const swim_step_t *start_step, DReferenceTrajectory *rt, vector<const DCoordinateSystem*> &wires, vector<DVector3> &shifts, vector<double> &errs, vector<double> &chisqv);
		fit_status_t LeastSquaresB(fit_type_t fit_type, DReferenceTrajectory *rt, double &chisq);
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

		int eventnumber;
		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
		double MAX_HIT_DIST;
		double CDC_Z_MIN;
		double CDC_Z_MAX;
		vector<DReferenceTrajectory*>rtv;
		DReferenceTrajectory *tmprt;
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
		TH1F *cdc_can_resi, *fdc_can_resi, *fdc_can_resi_cath;
};

#endif // _DTrack_factory_ALT2_

