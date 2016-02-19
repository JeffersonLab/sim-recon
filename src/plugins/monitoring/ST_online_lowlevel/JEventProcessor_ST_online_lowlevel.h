// $Id$
//
//    File: JEventProcessor_ST_online_lowlevel.h
// Created: Fri Jun 19 13:21:45 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_online_lowlevel_
#define _JEventProcessor_ST_online_lowlevel_

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
//***************** ROOT header files********************
//*******************************************************
#include <TMath.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
//***************** ST header files************************
//*********************************************************
#include "START_COUNTER/DSCDigiHit.h"
#include "START_COUNTER/DSCTDCDigiHit.h"
#include "START_COUNTER/DSCHit.h"
// ***************** DAQ header files (Triger Times)*********
//***********************************************************
#include "DAQ/Df250WindowRawData.h"
#include "DAQ/Df250PulseRawData.h"
#include "DAQ/Df250PulseTime.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250PulsePedestal.h"
#include "DAQ/Df250Config.h"
#include <TTAB/DTTabUtilities.h>
//****************************** Define some constants *********
//**************************************************************
const uint32_t NCHANNELS  = 30;     // number of scintillator paddles
const float_t  ADC_PT_RES = 0.0625; // fADC250 pulse time resolution (ns)
//const float_t  TDC_RES    = 0.0559; // f1TDC resolution (ns)
//***************** Declare Two Dimensional Histograms*************
//*****************************************************************
static TH2I *h2_st_adc_tdc_multi;
static TH2I *h2_st_adc_hit_multi;
static TH2I *h2_raw_pi_sector;
static TH2I *h2_raw_ped_sector;
static TH2I *h2_raw_pt_sector;
static TH2I *h2_adc_pp_sector;
static TH2I *h2_adc_pcpi_sector;
static TH2I *h2_adc_pt_sector;
static TH2I *h2_adc_ped_sector;
static TH2I *h2_st_time_vs_pcpi;
static TH2I *h2_st_time_vs_pp;
static TH2I *h2_raw_tdcTime_sec;
static TH2I *h2_tdcTime_sec;
static TH2I *h2_t_sec; 
static TH2I *h2_tTDC_sec; 
static TH2I *h2_tfADC_sec;
static TH2I *h2_dE_sec;
//***************** Declare One Dimensional Histograms****************
//********************************************************************
static TH1I *h1_adc_sec;
static TH1I *h1_tdc_sec;
static TH1I *h1_hit_sec;
static TH1I *st_num_events;
//***************** Declare Dynamic Arrays of waveform Histograms********
// &&&&&& define the Boolian variables needed &&&&&
//***********************************************************************
TH1I** h_amp_vs_sampl_chan = new TH1I*[NCHANNELS];
TH1I** h_amp_vs_sampl_chan150 = new TH1I*[NCHANNELS];
TH1I** h_amp_vs_sampl_chan1000 = new TH1I*[NCHANNELS];
TH1I** h_amp_vs_sampl_chan2000 = new TH1I*[NCHANNELS];
TH1I** h_amp_vs_sampl_chan3000 = new TH1I*[NCHANNELS];
TH1I** h_amp_vs_sampl_chan4000 = new TH1I*[NCHANNELS];
Bool_t bool_sec[NCHANNELS];
Bool_t bool_sec150[NCHANNELS];
Bool_t bool_sec1000[NCHANNELS];
Bool_t bool_sec2000[NCHANNELS];
Bool_t bool_sec3000[NCHANNELS];
Bool_t bool_sec4000[NCHANNELS];
Double_t adc_pp;
Double_t adc_t;
Double_t adc_ped;
Double_t adc_pi;
Double_t adc_pcpi;
Double_t tdc_t;
Double_t st_time;
//// Offsets from CCDB
vector<double> a_pedestals;
vector<double> adc_time_offsets;
vector<double> tdc_time_offsets;
double t_tdc_base;
double t_base;
double t_scale;
//**************Histograms limits Determination******************************
//***************************************************************************
const uint32_t ADC_MULTI_MIN  = 0.;      // Lower limit of adc multiplicity histogram
const uint32_t ADC_MULTI_MAX  = 31.;     // Upper limit of adc multiplicity histogram
const uint32_t ADC_MULTI_BINS = 31.;     // Number of bins in adc multiplicity histogram
const uint32_t TDC_MULTI_MIN  = 0.;      // Lower limit of tdc multiplicity histogram
const uint32_t TDC_MULTI_MAX  = 40.;     // Upper limit of tdc multiplicity histogram
const uint32_t TDC_MULTI_BINS = 40.;     // Number of bins in tdc multiplicity histogram
const uint32_t PI_MIN         = 0.;      // Lower limit of pulse integral histogram
const uint32_t PI_MAX         = 12000.;  // Upper limit of pulse integral histogram  
const uint32_t PI_BINS        = 240.;    // Number of bins in pulse integral histogram
const uint32_t PED_MIN        = 90.;     // Lower limit of pedestal histogram
const uint32_t PED_MAX        = 120.;    // Upper limit of pedestal histogram
const uint32_t PED_BINS       = 30.;     // Number of bins in pedestal histogram
const uint32_t PT_MIN         = 0.;      // Lower limit of pulse time histogram (ns)
const uint32_t PT_MAX         = 200.;    // Upper limit of pulse time histogram (ns)
const uint32_t PT_BINS        = 200.;    // Number of bins in pulse time histogram
const uint32_t TDC_DHIT_MIN   = 0.;      // Lower limit of tdc digihit time histogram 
const uint32_t TDC_DHIT_MAX   = 2200.;   // Upper limit of tdc digihit time histogram 
const uint32_t TDC_DHIT_BINS  = 220.;    // Number of bins in tdc digihit time histogram 
const float_t  T_HIT_MIN      =  -80.;  // Lower limit of hit time histogram (ns)
const float_t  T_HIT_MAX      =  80.;   // Upper limit of hit time histogram (ns)
const float_t  T_HIT_BINS     =  160.;   // Number of bins in hit time histogram 

class JEventProcessor_ST_online_lowlevel:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_online_lowlevel();
		~JEventProcessor_ST_online_lowlevel();
		const char* className(void){return "JEventProcessor_ST_online_lowlevel";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		
};

#endif // _JEventProcessor_ST_online_lowlevel_

