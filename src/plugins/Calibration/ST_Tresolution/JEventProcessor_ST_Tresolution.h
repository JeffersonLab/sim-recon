// $Id$
//
//    File: JEventProcessor_ST_Tresolution.h
// Created: Fri Jan  8 09:07:34 EST 2016
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_Tresolution_
#define _JEventProcessor_ST_Tresolution_

#include <JANA/JEventProcessor.h>
// ROOT header files
#include <TMath.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
// C++ header files
#include <stdint.h>
#include <vector>
#include <stdio.h>
// PID libraries
#include <PID/DEventRFBunch.h>
#include <PID/DParticleID.h>
#include <PID/DDetectorMatches.h>
#include <PID/DParticleID.h>
#include <PID/DChargedTrack.h>
// Tracking libraries
#include <TRACKING/DTrackTimeBased.h>
// RF libraries
#include <RF/DRFTDCDigiTime.h>
#include <RF/DRFTime_factory.h>
// ST libraries
#include <START_COUNTER/DSCHit.h>
#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <START_COUNTER/DSCTruthHit.h>
// TOF libraries
#include <TOF/DTOFHit.h>
// DAQ/EPICS libraries
#include <DAQ/DEPICSvalue.h>
// Translation table libraries
#include <TTAB/DTTabUtilities.h>
#include <TTAB/DTranslationTable.h>
//
const Int_t NCHANNELS = 30;
class JEventProcessor_ST_Tresolution:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_Tresolution();
		~JEventProcessor_ST_Tresolution();
		const char* className(void){return "JEventProcessor_ST_Tresolution";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		const DParticleID* dParticleID;
		DRFTime_factory *dRFTimeFactory;
		double z_target_center;  // Target center along z
		double dRFBunchPeriod;   // RF bunch period from CCDBr
		//////////////////////////////////
		double locSCHitTime,       locSCTrackFlightTime,   locFlightTimeCorrectedSCTime;
		double locTOFHitTime,      locTOFTrackFlightTime,  locFlightTimeCorrectedTOFTime;
		double locCenteredRFTime,  locCenterToVertexRFTime,  locVertexRFTime, locPeriod;
		double SC_RFShiftedTime;
		double locSCzIntersection;
		double locSCPropTime;
		double sc_pos_soss, sc_pos_eoss, sc_pos_eobs, sc_pos_eons;
		double  Corr_Time,Corr_Time_ss,Corr_Time_bs,Corr_Time_ns, Corr_Time_bn;
		double time_lower_limit, time_upper_limit, z_lower_limit, z_upper_limit;
		int NoBins_time, NoBins_z;
		double st_time, FlightTime, st_corr_FlightTime; 
		vector<vector<DVector3> >sc_pos;   // SC geometry vector
		vector<vector<DVector3> >sc_norm;  
		// Grab match track detector parameters
		DTOFHitMatchParams locTOFHitMatchParams;  // TOF
		DSCHitMatchParams  locSCHitMatchParams;   // SC
		// Declare event vertex vector
		DVector3 vertex;
		// Declare a vector which quantizes the point of the intersection of a charged particle 
		//   with a plane in the middle of the scintillator 
		DVector3 IntersectionPoint;
		// Declare a vector which quantizes the unit vector of the charged particle track traversing
	//   through the scintillator with its origin at the intersection point
		DVector3 IntersectionDir;
		// Grab the paramteres associated to a track matched to the ST
		vector<DSCHitMatchParams> st_params;
		// Declare r and z coord. of track vertex, sc array index
		double z_v, r_v;
		int sc_index;
		// Cuts
		double trackingFOMCut;
		double pim_pmag_cut;
		bool z_vertex_cut, r_vertex_cut;
		bool foundTOF, foundSC, foundSCandTOF;
		bool sc_match, sc_match_pid;
		// ******************* 2D histos **************
		TH2I **h2_CorrectedTime_z;
		//Define Calibration parameters variable called from CCDB
		vector<vector<double> >propagation_time_corr;


};

#endif // _JEventProcessor_ST_Tresolution_

