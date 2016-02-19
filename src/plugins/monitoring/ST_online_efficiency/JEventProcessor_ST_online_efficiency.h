// $Id$
//
//    File: JEventProcessor_ST_online_efficiency.h
// Created: Wed Jan 20 10:35:58 EST 2016
// Creator: mkamel (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_online_efficiency_
#define _JEventProcessor_ST_online_efficiency_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ST header files
#include "START_COUNTER/DSCHit.h"
#include "START_COUNTER/DSCDigiHit.h"
// PID libraries
#include "PID/DParticleID.h"
#include "PID/DChargedTrack.h"
#include <PID/DDetectorMatches.h>
#include <TRACKING/DTrackFitter.h>
//#include <PID/DDetectorMatches.h>
// Tracking libraries
#include "TRACKING/DTrackTimeBased.h"
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
// Define some constants
const double  RAD2DEG       = 57.29577951;      // Convert radians to degrees
const uint32_t  NCHANNELS     = 30;              // number of scintillator paddles
const uint32_t  Nof_ss_intervals = 5;
const uint32_t  Nof_bs_intervals = 3;
const uint32_t  Nof_ns_intervals = 4;
// Declare 1D tracking histos
static TH1D *h_N_trck_prd_All;
static TH1D *h_N_recd_hit_All;
static TH1D *h_N_trck_prd_ss;
static TH1D *h_N_recd_hit_ss;
static TH1D *h_N_trck_prd_bs;
static TH1D *h_N_recd_hit_bs;
static TH1D *h_N_trck_prd_ns;
static TH1D *h_N_recd_hit_ns;
static TH1D *h_ST_Eff_All;
static TH1D *h_ST_Eff_ss;
static TH1D *h_ST_Eff_bs;
static TH1D *h_ST_Eff_ns;


// Declare detection efficency counters
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

class JEventProcessor_ST_online_efficiency:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_online_efficiency();
		~JEventProcessor_ST_online_efficiency();
		const char* className(void){return "JEventProcessor_ST_online_efficiency";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		const DParticleID* dParticleID;
		DSCHitMatchParams  locSCHitMatchParams;   // SC
		vector<vector<DVector3> >sc_pos;   // SC geometry vector
		vector<vector<DVector3> >sc_norm;  
		double z_target_center;  // Target center along z
		double sc_pos_soss, sc_pos_eoss, sc_pos_eobs, sc_pos_eons;
		double ss_interval,bs_interval,ns_interval;
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
		double locSCzIntersection;
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
};

#endif // _JEventProcessor_ST_online_efficiency_

