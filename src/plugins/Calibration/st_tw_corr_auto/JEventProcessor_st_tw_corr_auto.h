// $Id$
//
//    File: JEventProcessor_st_tw_corr_auto.h
// Created: Mon Oct 26 10:35:45 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_st_tw_corr_auto_
#define _JEventProcessor_st_tw_corr_auto_

#include <JANA/JEventProcessor.h>
using namespace jana;
using namespace std;
using std::vector;
using std::string;
// ***************** C++ header files******************
//*****************************************************
#include <stdint.h>
#include <vector>
#include <stdio.h>
#include <fstream>
// ***************** ST header files*******************************
#include "START_COUNTER/DSCHit.h"
#include "START_COUNTER/DSCDigiHit.h"
// ***************** ROOT header files*****************************
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
// ******************* DAQ libraries******************************
#include <DAQ/Df250PulsePedestal.h>

class JEventProcessor_st_tw_corr_auto:public jana::JEventProcessor{
	public:
		JEventProcessor_st_tw_corr_auto();
		~JEventProcessor_st_tw_corr_auto();
		const char* className(void){return "JEventProcessor_st_tw_corr_auto";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

        
        // ***************** define constants*************************
        Int_t NCHANNELS;
        Double_t tdc_thresh_mV;
        Double_t tdc_gain_factor;
        Double_t adc_max_chan;
        Double_t adc_max_mV;
        Double_t adc_thresh_calc ;
        // ******************* 2D histos **************
        TH2I **h2_stt_vs_pp_chan;
	TH2I ** h2_st_corr_vs_pp;
        // ******************** 1D histos *******************
        TH1I **h_pp_chan;
        TH1I **h_stt_chan;
	TH1I **h1_st_corr_time;
	//Define Calibration parameters variable called from CCDB
	vector<vector<double> >timewalk_parameters;
	double USE_TIMEWALK_CORRECTION;
};

#endif // _JEventProcessor_st_tw_corr_auto_

