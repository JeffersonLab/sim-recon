// $Id$
//
//    File: JEventProcessor_ST_ZEff.h
// Created: Wed Aug 26 17:18:47 EDT 2015
// Creator: mkamel (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_ZEff_
#define _JEventProcessor_ST_ZEff_
// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
// ST header files
#include "START_COUNTER/DSCHit.h"
#include "START_COUNTER/DSCDigiHit.h"
// PID libraries
#include "PID/DParticleID.h"
#include "PID/DChargedTrack.h"
#include <PID/DDetectorMatches.h>
#include <TRACKING/DTrackFitter.h>
//#include <PID/DDetectorMatches.h>
// RF libraries
#include <RF/DRFTDCDigiTime.h>
#include <RF/DRFTime_factory.h>
// Tracking libraries
#include "TRACKING/DTrackTimeBased.h"
// TOF libraries
#include <TOF/DTOFHit.h>
#include <BCAL/DBCALHit.h>
// ROOT header files
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
// C++ header files
#include <stdint.h>
#include <vector>
#include <stdio.h>
#include "TObjArray.h"
#include "TMath.h"
//using std::vector;
//using std::string;
// Define some constants
const double  RAD2DEG       = 57.29577951;      // Convert radians to degrees
const uint32_t  NCHANNELS     = 30;              // number of scintillator paddles
const uint32_t  Nof_ss_intervals = 8;
const uint32_t  Nof_bs_intervals = 4;
const uint32_t  Nof_ns_intervals = 8;
// Declare 1D tracking histos
static TH1D *h_z_v;
static TH1D *h_fom;
static TH1D *h_N_Hit_in_track;
static TH1D *h1_qVectorSize;
static TH1D *h1_qVectorSize_ACuts;

static TH1D *h1_RFtime;
static TH1D *h1_Centered_RFtime;
static TH1D *h1_SC_ShiftedTime;
static TH1D *h1_tDiff;

static TH1D *h_N_trck_prd_z_ss[Nof_ss_intervals];
static TH1D *h_N_recd_hit_z_ss[Nof_ss_intervals];
static TH1D *h_N_recd_hit_z_ss_ACC[Nof_ss_intervals];

static TH1D *h_N_trck_prd_z_bs[Nof_bs_intervals];
static TH1D *h_N_recd_hit_z_bs[Nof_bs_intervals];
static TH1D *h_N_recd_hit_z_bs_ACC[Nof_bs_intervals];

static TH1D *h_N_trck_prd_z_ns[Nof_ns_intervals];
static TH1D *h_N_recd_hit_z_ns[Nof_ns_intervals];
static TH1D *h_N_recd_hit_z_ns_ACC[Nof_ns_intervals];

static TH1D *h_N_trck_prd_z_ss_eff[Nof_ss_intervals];
static TH1D *h_N_trck_prd_z_bs_eff[Nof_bs_intervals];
static TH1D *h_N_trck_prd_z_ns_eff[Nof_ns_intervals];
		static TH1D *h1_st_pred_id;
		// Declare 2D tracking histos
                
static TH2I *h2_z_vs_r;
static TH2I *h2_x_vs_y;






class JEventProcessor_ST_ZEff:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_ZEff();
		~JEventProcessor_ST_ZEff();
		const char* className(void){return "JEventProcessor_ST_ZEff";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		const DParticleID* dParticleID;
		double z_target_center;  // Target center along z
                 DRFTime_factory *dRFTimeFactory;
		vector<vector<DVector3> >sc_pos;   // SC geometry vector
		vector<vector<DVector3> >sc_norm;  
double Corr_Time_ns,Corr_Time_bs,Corr_Time_ss;
double SC_RFShiftedTime_ns,SC_RFShiftedTime_bs,SC_RFShiftedTime_ss;
double SC_RFShiftedTime;
double locVertexRFTime;


 double z_ss[Nof_ss_intervals],z_bs[Nof_bs_intervals],z_ns[Nof_ss_intervals]; 
 double incpt_ss,slope_ss,incpt_bs,slope_bs,incpt_ns,slope_ns,Bound1,Bound2; 
        /*
		// Grab match track detector parameters
		uint32_t N_trck_prd_All[NCHANNELS];
		uint32_t N_recd_hit_All[NCHANNELS];
		uint32_t N_trck_prd_ss[NCHANNELS];
		uint32_t N_recd_hit_ss[NCHANNELS];
		uint32_t N_trck_prd_bs[NCHANNELS];
		uint32_t N_recd_hit_bs[NCHANNELS];
		uint32_t N_trck_prd_ns[NCHANNELS];
		uint32_t N_recd_hit_ns[NCHANNELS];
		uint32_t N_recd_hit_All_nearby_plus[NCHANNELS];
		uint32_t N_recd_hit_All_nearby_minus[NCHANNELS];
		uint32_t N_recd_hit_All_nearby[NCHANNELS];
		uint32_t N_trck_prd_z_ss  [Nof_ss_intervals][NCHANNELS];
		uint32_t N_recd_hit_z_ss  [Nof_ss_intervals][NCHANNELS];
		uint32_t N_trck_prd_z_bs  [Nof_bs_intervals][NCHANNELS];
		uint32_t N_recd_hit_z_bs  [Nof_bs_intervals][NCHANNELS];
		uint32_t N_trck_prd_z_ns  [Nof_ns_intervals][NCHANNELS];
		uint32_t N_recd_hit_z_ns  [Nof_ns_intervals][NCHANNELS]; 
		DSCHitMatchParams  locSCHitMatchParams;   // SC
		double z_ss[Nof_ss_intervals],z_bs[Nof_bs_intervals],z_ns[Nof_ss_intervals]; 
		double theta_ss_max_left[Nof_ss_intervals];
		double theta_ss_max_right[Nof_ss_intervals];
		double theta_ss_min_left[Nof_ss_intervals];
		double theta_ss_min_right[Nof_ss_intervals];
		double smallest_ss_left[Nof_ss_intervals];
		double theta_ss_small[Nof_ss_intervals];
		double theta_ss_large[Nof_ss_intervals];
	
		double theta_bs_max_left[Nof_bs_intervals];
		double theta_bs_max_right[Nof_bs_intervals];
		double theta_bs_min_left[Nof_bs_intervals];
		double theta_bs_min_right[Nof_bs_intervals];
		double smallest_bs_left[Nof_bs_intervals];
		double theta_bs_small[Nof_bs_intervals];
		double theta_bs_large[Nof_bs_intervals];
	
		double theta_ns_max_left[Nof_ns_intervals];
		double theta_ns_max_right[Nof_ns_intervals];
		double theta_ns_min_left[Nof_ns_intervals];
		double theta_ns_min_right[Nof_ns_intervals];
		double smallest_ns_left[Nof_ns_intervals];
		double theta_ns_small[Nof_ns_intervals];
		double theta_ns_large[Nof_ns_intervals];
                bool theta_momentum_cut_ss[Nof_ss_intervals]; 
		bool theta_momentum_cut_bs[Nof_bs_intervals]; 	
		double locSCrIntersection;
		int sc_index;
		DVector3 locProjPos;
		// Declare a vector which quantizes the point of the intersection of a charged particle 
		//   with a plane in the middle of the scintillator 
		DVector3 IntersectionPoint;
		// Declare a vector which quantizes the unit vector of the charged particle track traversing
		//   through the scintillator with its origin at the intersection point
		DVector3 IntersectionDir;
		// Grab the paramteres associated to a track matched to the ST
		vector<DSCHitMatchParams> st_params;
		bool foundSC;
		bool sc_match, sc_match_pid;
		bool Barrel;
        */

		double sc_angle_corr;
		
        vector<vector<double> >propagation_time_corr;
	vector<vector<double> >PTC_Boundary;
};

#endif // _JEventProcessor_ST_ZEff_

