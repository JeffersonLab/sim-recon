// $Id$
//
//    File: JEventProcessor_ST_online_lowlevel.cc
// Created: Fri Jun 19 13:21:45 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//
#include "JEventProcessor_ST_online_lowlevel.h"
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
const uint32_t PI_MAX         = 40000.;  // Upper limit of pulse integral histogram  
const uint32_t PI_BINS        = 800.;    // Number of bins in pulse integral histogram
const uint32_t PED_MIN        = 90.;     // Lower limit of pedestal histogram
const uint32_t PED_MAX        = 130.;    // Upper limit of pedestal histogram
const uint32_t PED_BINS       = 40.;     // Number of bins in pedestal histogram
const uint32_t PT_MIN         = 0.;      // Lower limit of pulse time histogram (ns)
const uint32_t PT_MAX         = 200.;    // Upper limit of pulse time histogram (ns)
const uint32_t PT_BINS        = 200.;    // Number of bins in pulse time histogram
const uint32_t TDC_DHIT_MIN   = 0.;      // Lower limit of tdc digihit time histogram 
const uint32_t TDC_DHIT_MAX   = 4000.;   // Upper limit of tdc digihit time histogram 
const uint32_t TDC_DHIT_BINS  = 400.;    // Number of bins in tdc digihit time histogram 
const float_t  T_HIT_MIN      =  -200.;  // Lower limit of hit time histogram (ns)
const float_t  T_HIT_MAX      =  200.;   // Upper limit of hit time histogram (ns)
const float_t  T_HIT_BINS     =  400.;   // Number of bins in hit time histogram 
// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_online_lowlevel());
}
} // "C"


//------------------
// JEventProcessor_ST_online_lowlevel (Constructor)
//------------------
JEventProcessor_ST_online_lowlevel::JEventProcessor_ST_online_lowlevel()
{

}

//------------------
// ~JEventProcessor_ST_online_lowlevel (Destructor)
//------------------
JEventProcessor_ST_online_lowlevel::~JEventProcessor_ST_online_lowlevel()
{

}
bool DSCHit_fadc_cmp(const DSCDigiHit *a, const DSCDigiHit *b)
{
    if (a->sector == b->sector) return (a->pulse_time < b->pulse_time);
    return (a->sector < b->sector);
}

bool DSCHit_tdc_cmp(const DSCTDCDigiHit *a, const DSCTDCDigiHit *b)
{
    if (a->sector == b->sector) return (a->time < b->time);
    return (a->sector < b->sector);
}

bool DSCHit_thit_cmp(const DSCHit *a, const DSCHit *b)
{
    if (a->sector == b->sector) return (a->t < b->t);
    return (a->sector < b->sector);
}
//------------------
// init
//------------------
jerror_t JEventProcessor_ST_online_lowlevel::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
  
  //Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("st_lowlevel")->cd();

  st_num_events = new TH1I("st_num_events","ST Number of events",1, 0.5, 1.5);

  // 1D occupancy Histos
  h1_adc_sec = new TH1I("h1_adc_sec", "ST fADC250 DigiHit Occupancy; Channel Number; fADC250 Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);
  h1_tdc_sec = new TH1I("h1_tdc_sec", "ST TDC DigiHit Occupancy; Channel Number; TDC Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);
  h1_hit_sec = new TH1I("h1_hit_sec", "ST Hit Occupancy; Channel Number; Hit Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);

  // 2D Multiplicity Histos
  h2_st_adc_tdc_multi = new TH2I("h2_st_adc_tdc_multi", "ST Total Multiplicity: TDC vs ADC; f1TDC Multiplicity; fADC250 Multiplicity", TDC_MULTI_BINS, TDC_MULTI_MIN + 0.5, TDC_MULTI_MAX + 0.5, ADC_MULTI_BINS, ADC_MULTI_MIN + 0.5, ADC_MULTI_MAX + 0.5);
  h2_st_adc_hit_multi = new TH2I("h2_st_adc_hit_multi", "ST Total Multiplicity: HIT vs ADC; Hit Multiplicity; fADC250 Multiplicity", TDC_MULTI_BINS, TDC_MULTI_MIN + 0.5, TDC_MULTI_MAX + 0.5, ADC_MULTI_BINS, ADC_MULTI_MIN + 0.5, ADC_MULTI_MAX + 0.5);

  // 2D Raw ADC data + TDC time if there is a n adc hit Histos
  h2_raw_pi_sector  = new TH2I("h2_raw_pi_sector", "ST fADC250 Pulse Integral; Channel Number; Pulse Integral(au)", NCHANNELS, 0.5, NCHANNELS + 0.5, PI_BINS, PI_MIN, PI_MAX);
  h2_raw_ped_sector = new TH2I("h2_raw_ped_sector", "ST fADC250 Pulse Pedestal; Channel Number;Pulse Pedestal(au)", NCHANNELS, 0.5, NCHANNELS + 0.5, PED_BINS, PED_MIN, PED_MAX);
  h2_raw_pt_sector  = new TH2I("h2_raw_pt_sector", "ST fADC250 Pulse Time; Channel Number; Pulse Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, PT_BINS, PT_MIN, PT_MAX);
  h2_tdcTime_sec    = new TH2I("h2_tdcTime_sec", "ST TDC Sector vs Time (While ADC hit); Channel Number; TDC Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, TDC_DHIT_BINS, TDC_DHIT_MIN, TDC_DHIT_MAX);

  // 2D ADC/TDC data offsets applied 
  h2_adc_pp_sector   = new TH2I("h2_adc_pp_sector","Pulse peak vs sector;Channel number; Pulse Peak (channels) ",NCHANNELS, 0.5, NCHANNELS + 0.5,400,0,4000);
  h2_adc_pcpi_sector = new TH2I("h2_adc_pcpi_sector", "ST fADC250 Pedstal corrected Pulse Integral; Channel Number; fADC250 Pulse Integral (au)", NCHANNELS, 0.5, NCHANNELS + 0.5, PI_BINS, PI_MIN, PI_MAX);
  h2_adc_pt_sector   = new TH2I("h2_adc_pt_sector", "ST fADC250 Pulse Time; Channel Number; fADC250 Pulse Time (ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, 400, -100., 100);
  h2_adc_ped_sector  = new TH2I("h2_adc_ped_sector", "ST fADC250 Pedestal; Channel Number; fADC250 Pedestal", NCHANNELS, 0.5, NCHANNELS + 0.5, 100, 1000., 7000.);
  h2_st_time_vs_pcpi = new TH2I("h2_st_time_vs_pcpi","Ped Corrected Pulse Integral vs #delta (t_{TDC} - t_{ADC});Pulse Integral(cahnnels) ;#delta (t_{TDC} - t_{ADC}) (ns)",PI_BINS, PI_MIN, PI_MAX,80, -10., 10.);
  h2_st_time_vs_pp   = new TH2I("h2_st_time_vs_pp","Pulse Peak vs #delta (t_{TDC} - t_{ADC});Pulse Peak (cahnnels);#delta (t_{TDC} - t_{ADC}) (ns)",400,0,4000,80, -10., 10.);

  // Raw TDC data Histos
  h2_raw_tdcTime_sec = new TH2I("h2_raw_tdcTime_sec", "ST TDC Sector vs Time; Channel Number; TDC Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, TDC_DHIT_BINS, TDC_DHIT_MIN, TDC_DHIT_MAX);

  // Raw Hit data Histos
  h2_t_sec    = new TH2I("h2_t_sec", "ST Hit Sector vs Time (walk corrected); Channel Number; TDC Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, T_HIT_BINS, T_HIT_MIN, T_HIT_MAX); 
  h2_tTDC_sec = new TH2I("h2_tTDC_sec", "ST Hit Sector vs Time(No walk corrected); Channel Number; TDC Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, T_HIT_BINS, T_HIT_MIN, T_HIT_MAX);
  h2_tfADC_sec= new TH2I("h2_tfADC_sec", "ST Hit Sector vs ADC-Time; Channel Number; ADC Time(ns)", NCHANNELS, 0.5, NCHANNELS + 0.5, T_HIT_BINS, T_HIT_MIN, T_HIT_MAX); 
  h2_dE_sec   = new TH2I("h2_dE_sec", "ST Hit Sector vs dE; Channel Number; dE(GeV)", NCHANNELS, 0.5, NCHANNELS + 0.5, 100, 0.0, 0.01);
  //======================== Creat root folder for waveforms and cd to it================
  gDirectory->mkdir("waveforms")->cd();

  // Waveform Histos
  for(unsigned int i = 0; i < NCHANNELS; i++)
    {
      h_amp_vs_sampl_chan[i]     = new TH1I(Form("amp_vs_sampl_chan_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);
      h_amp_vs_sampl_chan150[i]     = new TH1I(Form("amp_vs_sampl_chan150_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);
      
      h_amp_vs_sampl_chan1000[i] = new TH1I(Form("amp_vs_sampl_chan1000_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);
      
      h_amp_vs_sampl_chan2000[i] = new TH1I(Form("amp_vs_sampl_chan2000_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);
      
      h_amp_vs_sampl_chan3000[i] = new TH1I(Form("amp_vs_sampl_chan3000_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);

      h_amp_vs_sampl_chan4000[i] = new TH1I(Form("amp_vs_sampl_chan4000_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Sample Number; fADC250 Pulse Height (au)", i+1, 0+12*i, 12+12*i), 100, 0, 100);
      // ========================Define Boolians Arrays to be false======================= 
      bool_sec[i] = false;
      bool_sec150[i] = false;
      bool_sec1000[i] = false;
      bool_sec2000[i] = false;
      bool_sec3000[i] = false;
      bool_sec4000[i] = false;
      // =================================================================================
    }
  // cd back to main directory
  main->cd();
  japp->RootUnLock();
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_online_lowlevel::brun(JEventLoop *eventLoop, int runnumber)
{
  // load scale factors
  map<string,double> scale_factors;
if (eventLoop->GetCalib("/START_COUNTER/digi_scales", scale_factors))
    jout << "Error loading /START_COUNTER/digi_scales !" << endl;
  // t_scale (SC_ADC_SCALE)
  if (scale_factors.find("SC_ADC_TSCALE") != scale_factors.end())
    t_scale = scale_factors["SC_ADC_TSCALE"];
  else
    jerr << "Unable to get SC_ADC_TSCALE from /START_COUNTER/digi_scales !" 
	 << endl;
  
  // load base time offset
  map<string,double> base_time_offset;
  // t_base (SC_BASE_TIME_OFFSET)
  if (eventLoop->GetCalib("/START_COUNTER/base_time_offset",base_time_offset))
    jout << "Error loading /START_COUNTER/base_time_offset !" << endl;
  if (base_time_offset.find("SC_BASE_TIME_OFFSET") != base_time_offset.end())
    t_base = base_time_offset["SC_BASE_TIME_OFFSET"];
  else
    jerr << "Unable to get SC_BASE_TIME_OFFSET from /START_COUNTER/base_time_offset !" << endl;
  // t_tdc_base (SC_TDC_BASE_TIME_OFFSET)
  if (base_time_offset.find("SC_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
    t_tdc_base = base_time_offset["SC_TDC_BASE_TIME_OFFSET"];
  else
    jerr << "Unable to get SC_BASE_TIME_OFFSET from /START_COUNTER/base_time_offset !" << endl;
  // load constant tables
  // a_pedestals (pedestals)
  if (eventLoop->GetCalib("/START_COUNTER/pedestals", a_pedestals))
    jout << "Error loading /START_COUNTER/pedestals !" << endl;
  // adc_time_offsets (adc_timing_offsets)
  if (eventLoop->GetCalib("/START_COUNTER/adc_timing_offsets", adc_time_offsets))
    jout << "Error loading /START_COUNTER/adc_timing_offsets !" << endl;
  // tdc_time_offsets (tdc_timing_offsets)
  if (eventLoop->GetCalib("/START_COUNTER/tdc_timing_offsets", tdc_time_offsets)) jout << "Error loading /START_COUNTER/tdc_timing_offsets !" << endl;
  // This is called whenever the run number changes
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_online_lowlevel::evnt(JEventLoop *loop, int eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

  // Get the data objects first so we minimize the time we hold the ROOT mutex lock
  vector<const DSCDigiHit*> dscdigihits;                // ST fADC250 DigiHits
  vector<const DSCTDCDigiHit*> dsctdcdigihits;          // ST f1TDC DigiHits
  vector<const DSCHit*> dschits;                        // ST hits
  const DTTabUtilities*           TTabUtilities = NULL;
  loop->Get(dscdigihits);
  loop->Get(dsctdcdigihits);
  loop->Get(dschits);
  loop->GetSingle(TTabUtilities);
  uint32_t ADC_hits       = dscdigihits.size();
  uint32_t TDC_hits       = dsctdcdigihits.size();
  uint32_t Hits           = dschits.size();

  // Lock ROOT mutex so other threads won't interfere 
  japp->RootWriteLock();
  if( (dscdigihits.size()>0) || (dsctdcdigihits.size()>0) || (dschits.size()>0) )
    st_num_events->Fill(1);

  //Fill 2D multiplicity histos
  h2_st_adc_tdc_multi->Fill(TDC_hits, ADC_hits);
  h2_st_adc_hit_multi->Fill(Hits, ADC_hits);
//========================== DSCDigiHits ** ADC Hits=================
//*******************************************************************
  for(uint32_t i = 0; i < ADC_hits; i++)
     {
       // Sort the ADC hits
       //   sort(dscdigihits.begin(), dscdigihits.end(), DSCHit_fadc_cmp);
       // const DSCDigiHit *adc_dhit = dscdigihits[i];
//=========== Wave Form online O'Scope =====================
// Define some objects
       
       const Df250PulseIntegral *pulseintegral = NULL;
       const Df250PulsePedestal *pulsepedestal = NULL;
       const Df250WindowRawData *windowrawdata = NULL;

      // Define a vector to store the adc samples
       vector<uint16_t> samples;  // Declare empty vector
      // Get the sector values
      Int_t hit_sector_adc       = dscdigihits[i]->sector;
      Int_t hit_sector_adc_index = hit_sector_adc - 1;
      
	  dscdigihits[i]->GetSingle(pulseintegral);
	  dscdigihits[i]->GetSingle(pulsepedestal);
	  dscdigihits[i]->GetSingle(windowrawdata);
	  // Obtain the pedestal and raw window data
	  if (pulseintegral) 
	    {
	      pulseintegral->GetSingle(pulsepedestal);
	      pulseintegral->GetSingle(windowrawdata);
	    }
	  // Histogram the raw window data
	  if (windowrawdata)
	    {  const Df250PulsePedestal* PPobj = NULL;
	      dscdigihits[i]->GetSingle(PPobj);
	      if (PPobj != NULL)
		{
		  if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
		}
	      adc_pp   = PPobj->pulse_peak;
	      if ((100 < adc_pp) && (adc_pp <= 150))
		{
		  if (!bool_sec150[hit_sector_adc_index])
		    {
		      bool_sec150[hit_sector_adc_index] = true;
		      if (bool_sec150[hit_sector_adc_index])
		      {
			
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    h_amp_vs_sampl_chan150[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}
	      
	      if ((150 < adc_pp) && (adc_pp <= 1000))
		{
		  if (!bool_sec[hit_sector_adc_index])
		    {
		      bool_sec[hit_sector_adc_index] = true;
		      if (bool_sec[hit_sector_adc_index])
		      {
			
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    h_amp_vs_sampl_chan[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}
	      
	      if ( (1000 < adc_pp) && (adc_pp <= 2000))
		{
		  if (!bool_sec1000[hit_sector_adc_index])
		    {
		      bool_sec1000[hit_sector_adc_index] = true;
		      if (bool_sec1000[hit_sector_adc_index])
		      {
			
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    
			    h_amp_vs_sampl_chan1000[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}

	      if ( (2000 < adc_pp) && (adc_pp <= 3000))
		{
		  if (!bool_sec2000[hit_sector_adc_index])
		    {
		      bool_sec2000[hit_sector_adc_index] = true;
		      if (bool_sec2000[hit_sector_adc_index])
		      {
			
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    
			    h_amp_vs_sampl_chan2000[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}

	      if ( (3000 < adc_pp) && (adc_pp <= 4000))
		{
		  if (!bool_sec3000[hit_sector_adc_index])
		    {
		      bool_sec3000[hit_sector_adc_index] = true;
		      if (bool_sec3000[hit_sector_adc_index])
		      {
			
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    
			    h_amp_vs_sampl_chan3000[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}

	      if ( (4000 < adc_pp))
		{
		  if (!bool_sec4000[hit_sector_adc_index])
		    {
		      bool_sec4000[hit_sector_adc_index] = true;
		      if (bool_sec4000[hit_sector_adc_index])
		      {
			for (uint32_t j = 0; j < windowrawdata->samples.size(); j++)
			  {
			    samples.push_back(windowrawdata->samples[j]);
			    h_amp_vs_sampl_chan4000[hit_sector_adc_index]->Fill(j, samples[j]);
			  }
		      }
		    }
		}
	    } // Windowrawdata cut
	  // ================= Get the raw data variables from ADC digihit object==========
	  int hit_channel       = dscdigihits[i]->sector - 1;             // channel hit (ranging from 0 to NCHANNELS-1)
	  int adc_sector         = dscdigihits[i]->sector ;                // channel hit (ranging from 1 to NCHANNELS)   
	  uint32_t avg_pedestal  = dscdigihits[i]->pedestal;               // average pedestal (should be around 100 chan1)
	  uint32_t pulse_time    = dscdigihits[i]->pulse_time*ADC_PT_RES;  // converted pulse time to ns
	  uint32_t pulse_integral= dscdigihits[i]->pulse_integral;         // pulse integral
	  //Occupancy Histo
	  h1_adc_sec->Fill(adc_sector);
	  // Fill the 2D histo
	  h2_raw_pi_sector->Fill(adc_sector,pulse_integral);
	  h2_raw_ped_sector->Fill(adc_sector,avg_pedestal);
	  h2_raw_pt_sector->Fill(adc_sector,pulse_time);
	  //Apply the calibration constants=================================================
	  const Df250PulsePedestal* PPobj = NULL;
	  dscdigihits[i]->GetSingle(PPobj);
	  if (PPobj != NULL)
	    {
	      if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
	    }
       
       // Initialize pedestal to one found in CCDB, but override it
       // with one found in event if is available
       double pedestal = a_pedestals[hit_channel];
       const Df250PulseIntegral *PIobj = NULL;
       dscdigihits[i]->GetSingle(PIobj);
       if (PIobj != NULL) 
	 {
	   double single_sample_ped = (double)PIobj->pedestal;
	   double nsamples_integral = (double)PIobj->nsamples_integral;
	   double nsamples_pedestal = (double)PIobj->nsamples_pedestal;
	   pedestal = single_sample_ped * nsamples_integral/nsamples_pedestal;
	 }      	
      
       // Apply calibration constants here
       adc_ped  = pedestal;
       adc_pi   =  dscdigihits[i]->pulse_integral;
       adc_pcpi = adc_pi - adc_ped;
       if(PPobj != NULL)
	 adc_pp   = PPobj->pulse_peak;
       else
	 adc_pp =0;
       adc_t    =  dscdigihits[i]->pulse_time * t_scale  - adc_time_offsets[hit_channel] + t_base; // Convert to ns
       //cout << "adc_t =  "<< adc_t << endl;
       //Fill 2D Histos
       h2_adc_pp_sector->Fill(adc_sector,adc_pp);
       h2_adc_pcpi_sector->Fill(adc_sector,adc_pcpi);
       h2_adc_pt_sector->Fill(adc_sector,adc_t);
       h2_adc_ped_sector->Fill(adc_sector,adc_ped);
	  // Aquire the TDC DigiHits******* get the tdc hits when there is an adc hit*****
       	  for(uint32_t i = 0; i < TDC_hits; i++)
	    {
	      //sort(dsctdcdigihits.begin(), dsctdcdigihits.end(), DSCHit_tdc_cmp);
	      const DSCTDCDigiHit *tdc_dhit = dsctdcdigihits[i];
	      float tdc_dhit_time = TTabUtilities->Convert_DigiTimeToNs_F1TDC(tdc_dhit);//tdc_dhit->time*TDC_RES;
	      int tdc_sector = tdc_dhit->sector;
	      if (adc_sector == tdc_sector)
		{
		  h2_tdcTime_sec->Fill(tdc_sector,tdc_dhit_time);

		  tdc_t = TTabUtilities->Convert_DigiTimeToNs_F1TDC(tdc_dhit) - tdc_time_offsets[tdc_sector] + t_tdc_base;
		  // cout << "tdc_t =  "<< tdc_t << endl;

		  st_time = tdc_t - adc_t;
		  //cout << "st_time =  "<< st_time  << endl;

		  h2_st_time_vs_pcpi->Fill(adc_pcpi,st_time);
		  h2_st_time_vs_pp->Fill(adc_pp,st_time);

		}
	    }// End TDC loop
     }// End ADC loop
//========================== DSCTDCDigiHits ** TDC Hits=================
//*******************************************************************
  for(uint32_t i = 0; i < TDC_hits; i++)
    {
      const DSCTDCDigiHit *tdc_dhit = dsctdcdigihits[i];
      //float tdc_dhit_time = tdc_dhit->time*TDC_RES;
      float tdc_dhit_time = TTabUtilities->Convert_DigiTimeToNs_F1TDC(tdc_dhit);
      int tdc_sec = tdc_dhit->sector;
      //Fill tdc occupancy histo
      h1_tdc_sec->Fill(tdc_sec);
      //Fill 2D histo
      h2_raw_tdcTime_sec->Fill(tdc_sec,tdc_dhit_time);

    }// End TDC loop

//==========================  DSCHits  ** Hits after hit factory=================
//*******************************************************************
  for(uint32_t i = 0; i < dschits.size(); i++) 
    {
      const DSCHit *hit = dschits[i];
      int hit_sector    = hit->sector;
      float dE          = hit->dE;          // Energy loss in GeV
      float t           = hit->t;           // best time (walk-corrected tdc)
      float t_TDC       = hit->t_TDC;       // time from TDC, no walk correction
      float t_fADC      = hit->t_fADC;      // time from fADC
      
      h1_hit_sec->Fill(hit_sector);
      //Fill 2D histos
      h2_t_sec->Fill(hit_sector,t);         // walk corrected tdc time vs sector
      h2_tTDC_sec->Fill(hit_sector,t_TDC);  // No walk correction TDC time vs sector
      h2_tfADC_sec->Fill(hit_sector,t_fADC);// ADC time vs Sector
      h2_dE_sec->Fill(hit_sector,dE);       // dE vs sector 
      
    }// End Hit loop

  // Lock ROOT mutex so other threads won't interfere 
  japp->RootUnLock();
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_online_lowlevel::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_online_lowlevel::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}
