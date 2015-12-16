// $Id$
//
//    File: JEventProcessor_PS_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>

#include "JEventProcessor_PS_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include <PAIR_SPECTROMETER/DPSDigiHit.h>
#include <PAIR_SPECTROMETER/DPSHit.h>
#include <PAIR_SPECTROMETER/DPSGeometry.h>
#include <DAQ/Df250PulsePedestal.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

// number of PS columns (tiles) per arm
const int Ncols = DPSGeometry::NUM_FINE_COLUMNS; // 145
const int Narms = DPSGeometry::NUM_ARMS; // 2
//
const int NmultBins = 10; //number of bins for multiplicity histograms
//
double counts_cut = 500.0;
// root hist pointers
static TH1I *ps_num_events;
static TH1I *hHit_NHits;
static TH1I *hHit_Arm;
static TH2I *hHit_NHitsVsArm;
static TH1I *hHit_Occupancy[Narms];
static TH1I *hHit_Energy[Narms];
static TH2I *hHit_EnergyVsColumn[Narms];
static TH1I *hHit_Integral[Narms];
static TH2I *hHit_IntegralVsColumn[Narms];
//static TH1I *hHit_Npix[Narms];
static TH1I *hHit_Time[Narms];
static TH2I *hHit_TimeVsColumn[Narms];
//
static TH1I *hDigiHit_NHits;
static TH1I *hDigiHit_Arm;
static TH2I *hDigiHit_NHitsVsArm;
static TH1I *hDigiHit_NSamplesPedestal[Narms];
static TH1I *hDigiHit_Pedestal[Narms];
static TProfile *hDigiHit_PedestalVsColumn[Narms];
static TH1I *hDigiHit_QualityFactor[Narms];
static TH1I *hDigiHit_Occupancy[Narms];
static TH1I *hDigiHit_RawPeak[Narms];
static TH2I *hDigiHit_RawPeakVsColumn[Narms];
static TH1I *hDigiHit_RawIntegral[Narms];
static TH2I *hDigiHit_RawIntegralVsColumn[Narms];
static TH1I *hDigiHit_NSamplesIntegral[Narms];
static TH2I *hDigiHit_PeakVsColumn[Narms];
static TH2I *hDigiHit_IntegralVsPeak[Narms];
static TH2I *hDigiHit_IntegralVsColumn[Narms];
static TH1I *hDigiHit_PulseTime[Narms];
static TH1I *hDigiHit_Time[Narms];
static TH2I *hDigiHit_TimeVsColumn[Narms];
static TH2I *hDigiHit_TimeVsIntegral[Narms];
//
static TH1I *hDigiHit_NHits_cut;
static TH1I *hDigiHit_Arm_cut;
static TH2I *hDigiHit_NHitsVsArm_cut;
static TH1I *hDigiHit_Occupancy_cut[Narms];
static TH1I *hDigiHit_Time_cut[Narms];
static TH2I *hDigiHit_TimeVsColumn_cut[Narms];
static TH2I *hDigiHit_TimeVsQF_cut[Narms];
//----------------------------------------------------------------------------------

// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_PS_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_PS_online::JEventProcessor_PS_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_PS_online::~JEventProcessor_PS_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_PS_online::init(void) {

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  // create root folder for ps and cd to it, store main dir
  TDirectory *mainDir = gDirectory;
  TDirectory *psDir = gDirectory->mkdir("PS");
  psDir->cd();

  // book hists
  ps_num_events = new TH1I("ps_num_events","PS Number of events",1,0.5,1.5);
  // fADC250 hit-level hists (after calibration)	  
  TDirectory *hitDir = gDirectory->mkdir("Hit"); hitDir->cd();
  hHit_NHits = new TH1I("Hit_NHits","PS hit multiplicity;hits;events",NmultBins,0.5,0.5+NmultBins);
  hHit_Arm = new TH1I("Hit_Arm","PS arm;arm;hits",Narms,-0.5,-0.5+Narms);
  hHit_NHitsVsArm    = new TH2I("Hit_NHitsVsArm","PS hit multiplicity vs. arm;arm;hits",Narms,-0.5,-0.5+Narms,NmultBins,0.5,0.5+NmultBins);
  TString arm_str[] = {"Left","Right"};
  for (int i=0;i<Narms;i++) {
    gDirectory->mkdir(arm_str[i]+"Arm")->cd();
    TString strN = "_" + arm_str[i] + "Arm";
    TString strT = ", " + arm_str[i] + " arm";
    hHit_Occupancy[i]    = new TH1I("Hit_Occupancy"+strN,"PS occupancy"+strT+";column (tile);hits / column",Ncols,0.5,0.5+Ncols);
    hHit_Energy[i] = new TH1I("Hit_Energy"+strN,"PS energy"+strT+";energy [GeV];hits / column",120,0.5,6.5);  
    hHit_EnergyVsColumn[i] = new TH2I("Hit_EnergyVsColumn"+strN,"PS energy vs. column"+strT+";column (tile);energy [GeV]",Ncols,0.5,0.5+Ncols,120,0.5,6.5);
    hHit_Integral[i]    = new TH1I("Hit_Integral"+strN,"PS fADC pulse integral"+strT+";pulse integral;hits",1000,0.0,30000.0);
    hHit_IntegralVsColumn[i]    = new TH2I("Hit_IntegralVsColumn"+strN,"PS fADC pulse integral vs. column"+strT+";column (tile);pulse integral",Ncols,0.5,0.5+Ncols,1000,0.0,30000.0);
    //hHit_Npix[i]    = new TH1I("Hit_Npix"+strN,"PS fADC number of pixels"+strT+";pixels;hits",200,0,100);
    hHit_Time[i]    = new TH1I("Hit_Time"+strN,"PS fADC time"+strT+";time [ns];hits / 400 ps",1000,-200.0,200.0);
    hHit_TimeVsColumn[i]    = new TH2I("Hit_TimeVsColumn"+strN,"PS fADC time vs. column"+strT+";column (tile);time [ns]",Ncols,0.5,0.5+Ncols,1000,-200.0,200.0);
    hitDir->cd();
  }
  // fADC250 digihit-level hists
  psDir->cd();
  TDirectory *digihitDir = gDirectory->mkdir("DigiHit"); digihitDir->cd();
  hDigiHit_NHits = new TH1I("DigiHit_NHits","PS fADC hit multiplicity;raw hits;events",NmultBins,0.5,0.5+NmultBins);
  hDigiHit_Arm = new TH1I("DigiHit_Arm","PS arm;arm;raw hits",Narms,-0.5,-0.5+Narms);
  hDigiHit_NHitsVsArm    = new TH2I("DigiHit_NHitsVsArm","PS hit multiplicity vs. arm;arm;raw hits",Narms,-0.5,-0.5+Narms,NmultBins,0.5,0.5+NmultBins);
  hDigiHit_NHits_cut = new TH1I("DigiHit_NHits_cut","PS fADC hit multiplicity (> 500 ADC integral counts);raw hits;events",NmultBins,0.5,0.5+NmultBins);
  hDigiHit_Arm_cut = new TH1I("DigiHit_Arm_cut","PS arm (> 500 ADC integral counts);arm;raw hits",Narms,-0.5,-0.5+Narms);
  hDigiHit_NHitsVsArm_cut    = new TH2I("DigiHit_NHitsVsArm_cut","PS hit multiplicity vs. arm (> 500 ADC integral counts);arm;raw hits",Narms,-0.5,-0.5+Narms,NmultBins,0.5,0.5+NmultBins);
  for (int i=0;i<Narms;i++) {
    gDirectory->mkdir(arm_str[i]+"Arm")->cd();
    TString strN = "_" + arm_str[i] + "Arm";
    TString strT = ", " + arm_str[i] + " arm";
    hDigiHit_NSamplesPedestal[i]    = new TH1I("DigiHit_NSamplesPedestal"+strN,"PS fADC pedestal samples"+strT+";pedestal samples;raw hits",50,-0.5,49.5);
    hDigiHit_Pedestal[i] = new TH1I("DigiHit_Pedestal"+strN,"PS fADC pedestals"+strT+";pedestal [fADC counts];raw hits",200,0.0,200.0);
    hDigiHit_PedestalVsColumn[i] = new TProfile("DigiHit_PedestalVsColumn"+strN,"PS pedestal vs. column"+strT+";column (tile);average pedestal [fADC counts]",Ncols,0.5,0.5+Ncols,"s"); 
    hDigiHit_QualityFactor[i] = new TH1I("DigiHit_QualityFactor"+strN,"PS fADC quality factor"+strT+";quality factor;raw hits",4,-0.5,3.5);
    hDigiHit_Occupancy[i]    = new TH1I("DigiHit_Occupancy"+strN,"PS fADC hit occupancy"+strT+";column (tile);raw hits / column",Ncols,0.5,0.5+Ncols);
    hDigiHit_RawPeak[i]    = new TH1I("DigiHit_RawPeak"+strN,"PS fADC pulse peak (raw)"+strT+";pulse peak (raw);raw hits",410,0.0,4100.0);
    hDigiHit_RawPeakVsColumn[i]    = new TH2I("DigiHit_RawPeakVsColumn"+strN,"PS fADC pulse peak (raw) vs. column"+strT+";column (tile);pulse peak (raw)",Ncols,0.5,0.5+Ncols,410,0.0,4100.0);
    hDigiHit_RawIntegral[i]    = new TH1I("DigiHit_RawIntegral"+strN,"PS fADC pulse integral (raw)"+strT+";pulse integral (raw);raw hits",1000,0.0,30000.0);
    hDigiHit_RawIntegralVsColumn[i]    = new TH2I("DigiHit_RawIntegralVsColumn"+strN,"PS fADC pulse integral (raw) vs. column"+strT+";column (tile);pulse integral (raw)",Ncols,0.5,0.5+Ncols,1000,0.0,30000.0);
    hDigiHit_NSamplesIntegral[i]    = new TH1I("DigiHit_NSamplesIntegral"+strN,"PS fADC integral samples"+strT+";integral samples;raw hits",60,-0.5,59.5);
    hDigiHit_PeakVsColumn[i]    = new TH2I("DigiHit_PeakVsColumn"+strN,"PS fADC pulse peak vs. column"+strT+";column (tile);pulse peak",Ncols,0.5,0.5+Ncols,410,0.0,4100.0);
    hDigiHit_IntegralVsPeak[i]    = new TH2I("DigiHit_IntegralVsPeak"+strN,"PS fADC pulse integral vs. peak"+strT+";pulse peak;pulse integral",410,0.0,4100.0,1000,0.0,30000.0);
    hDigiHit_IntegralVsColumn[i]    = new TH2I("DigiHit_IntegralVsColumn"+strN,"PS fADC pulse integral vs. column"+strT+";column (tile);pulse integral",Ncols,0.5,0.5+Ncols,1000,0.0,30000.0);
    hDigiHit_PulseTime[i]    = new TH1I("DigiHit_PulseTime"+strN,"PS fADC pulse time"+strT+";pulse time [62.5 ps];raw hits",2000,0.0,6500.0);
    hDigiHit_Time[i]    = new TH1I("DigiHit_Time"+strN,"PS fADC pulse time"+strT+";pulse time [ns];raw hits / 400 ps",1000,0.0,400.0);
    hDigiHit_TimeVsColumn[i]    = new TH2I("DigiHit_TimeVsColumn"+strN,"PS fADC pulse time vs. column"+strT+";column (tile);pulse time [ns]",Ncols,0.5,0.5+Ncols,1000,0.0,400.0);
    hDigiHit_TimeVsIntegral[i]    = new TH2I("DigiHit_TimeVsIntegral"+strN,"PS fADC pulse time vs. integral"+strT+";pulse integral;pulse time [ns]",1000,0.0,30000.0,1000,0.0,400.0);
    // digihit-level hists after cut on ADC integral counts
    hDigiHit_Occupancy_cut[i]    = new TH1I("DigiHit_Occupancy_cut"+strN,"PS fADC hit occupancy (> 500 ADC integral counts)"+strT+";column (tile);raw hits / column",Ncols,0.5,0.5+Ncols);
    hDigiHit_Time_cut[i]    = new TH1I("DigiHit_Time_cut"+strN,"PS fADC pulse time (> 500 ADC integral counts)"+strT+";pulse time [ns];raw hits / 400 ps",1000,0.0,400.0);
    hDigiHit_TimeVsColumn_cut[i]    = new TH2I("DigiHit_TimeVsColumn_cut"+strN,"PS fADC pulse time vs. column (> 500 ADC integral counts)"+strT+";column (tile);pulse time [ns]",Ncols,0.5,0.5+Ncols,1000,0.0,400.0);
    hDigiHit_TimeVsQF_cut[i]    = new TH2I("DigiHit_TimeVsQF_cut"+strN,"PS fADC pulse time vs. quality factor (> 500 ADC integral counts)"+strT+";fADC quality factor;pulse time [ns]",4,-0.5,3.5,1000,0.0,400.0);
    digihitDir->cd();
  }
  // back to main dir
  mainDir->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_PS_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  // extract the PS geometry
  vector<const DPSGeometry*> psGeomVect;
  eventLoop->Get(psGeomVect);
  if (psGeomVect.size() < 1)
    return OBJECT_NOT_AVAILABLE;
  const DPSGeometry& psGeom = *(psGeomVect[0]);
  // get photon energy bin lows for variable-width energy binning
  double Elows[Narms][Ncols+1];
  double cols[Ncols+1];
  for (int i=0;i<Ncols;i++) {
    Elows[0][i] = psGeom.getElow(0,i+1); 
    Elows[1][i] = psGeom.getElow(1,i+1); 
    cols[i] = 0.5+i;
  }
  // add the upper limits
  Elows[0][Ncols] = psGeom.getEhigh(0,Ncols); 
  Elows[1][Ncols] = psGeom.getEhigh(1,Ncols);
  cols[Ncols] = 0.5+Ncols;
  //  
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
  // set variable-width energy bins if histogram is empty
  for (int i=0;i<Narms;i++) {
    if (hHit_Energy[i]->GetEntries()==0) hHit_Energy[i]->SetBins(Ncols,Elows[i]);
    if (hHit_EnergyVsColumn[i]->GetEntries()==0) hHit_EnergyVsColumn[i]->SetBins(Ncols,cols,Ncols,Elows[i]);
  }
  japp->RootUnLock(); //RELEASE ROOT LOCK!!
  //
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_PS_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  // get data for PS
  vector<const DPSDigiHit*> digihits;
  eventLoop->Get(digihits); 
  vector<const DPSHit*> hits;
  eventLoop->Get(hits);
  // cache pulse pedestal objects
  map< const DPSDigiHit*, const Df250PulsePedestal* > pp_cache;
  for(unsigned int i=0; i < digihits.size(); i++) {
    const Df250PulsePedestal* pulsePed = NULL;
    digihits[i]->GetSingle(pulsePed);
    pp_cache[digihits[i]] = pulsePed;
  }

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  if(digihits.size()>0) ps_num_events->Fill(1);
  
  // fill digihit hists
  int NDigiHits[] = {0,0}; 
  int NDigiHits_cut[] = {0,0}; 
  hDigiHit_NHits->Fill(digihits.size());
  for(unsigned int i=0; i < digihits.size(); i++) {
    int arm = digihits[i]->arm;
    double ped = digihits[i]->pedestal/digihits[i]->nsamples_pedestal;
    hDigiHit_NSamplesPedestal[arm]->Fill(digihits[i]->nsamples_pedestal);
    hDigiHit_Pedestal[arm]->Fill(ped);
    const Df250PulsePedestal* pulsePed = pp_cache[digihits[i]];
    double peak = -999.0;
    if (pulsePed) peak = pulsePed->pulse_peak; 
    hDigiHit_RawPeak[arm]->Fill(peak);
    if (ped==0.0||peak==0.0) continue;
    hDigiHit_PedestalVsColumn[arm]->Fill(digihits[i]->column,ped);
    NDigiHits[arm]++;
    hDigiHit_Arm->Fill(arm);
    hDigiHit_Occupancy[arm]->Fill(digihits[i]->column);
    hDigiHit_RawPeakVsColumn[arm]->Fill(digihits[i]->column,peak);
    hDigiHit_RawIntegral[arm]->Fill(digihits[i]->pulse_integral);
    hDigiHit_RawIntegralVsColumn[arm]->Fill(digihits[i]->column,digihits[i]->pulse_integral);
    hDigiHit_NSamplesIntegral[arm]->Fill(digihits[i]->nsamples_integral);
    hDigiHit_PeakVsColumn[arm]->Fill(digihits[i]->column,peak-ped);
    double PI = digihits[i]->pulse_integral-digihits[i]->nsamples_integral*ped; // pedestal-subtracted pulse integral
    hDigiHit_IntegralVsColumn[arm]->Fill(digihits[i]->column,PI);
    hDigiHit_IntegralVsPeak[arm]->Fill(peak-ped,PI);
    hDigiHit_PulseTime[arm]->Fill(digihits[i]->pulse_time);
    double t_ns = 0.0625*digihits[i]->pulse_time;
    hDigiHit_Time[arm]->Fill(t_ns);
    hDigiHit_TimeVsColumn[arm]->Fill(digihits[i]->column,t_ns);
    hDigiHit_TimeVsIntegral[arm]->Fill(PI,t_ns);
    hDigiHit_QualityFactor[arm]->Fill(digihits[i]->QF);
    if (PI>counts_cut)  { 
      NDigiHits_cut[arm]++;
      hDigiHit_Arm_cut->Fill(arm);
      hDigiHit_Occupancy_cut[arm]->Fill(digihits[i]->column);
      hDigiHit_Time_cut[arm]->Fill(t_ns);
      hDigiHit_TimeVsColumn_cut[arm]->Fill(digihits[i]->column,t_ns);
      hDigiHit_TimeVsQF_cut[arm]->Fill(digihits[i]->QF,t_ns);
    }
  }  
  hDigiHit_NHitsVsArm->Fill(0.,NDigiHits[0]); hDigiHit_NHitsVsArm->Fill(1.,NDigiHits[1]);
  hDigiHit_NHits_cut->Fill(NDigiHits_cut[0]+NDigiHits_cut[1]);
  hDigiHit_NHitsVsArm_cut->Fill(0.,NDigiHits_cut[0]); hDigiHit_NHitsVsArm_cut->Fill(1.,NDigiHits_cut[1]);
  // fill hit hists
  int NHits[] = {0,0}; 
  hHit_NHits->Fill(hits.size());
  for(unsigned int i=0; i<hits.size(); i++) {
    int arm = hits[i]->arm;
    NHits[arm]++;
    hHit_Arm->Fill(arm);
    hHit_Occupancy[arm]->Fill(hits[i]->column);
    hHit_Energy[arm]->Fill(hits[i]->E);
    hHit_EnergyVsColumn[arm]->Fill(hits[i]->column,hits[i]->E);
    hHit_Integral[arm]->Fill(hits[i]->integral);
    hHit_IntegralVsColumn[arm]->Fill(hits[i]->column,hits[i]->integral);
    //hHit_Npix[arm]->Fill(hits[i]->npix_fadc);
    hHit_Time[arm]->Fill(hits[i]->t);
    hHit_TimeVsColumn[arm]->Fill(hits[i]->column,hits[i]->t);
  }
  hHit_NHitsVsArm->Fill(0.,NHits[0]); hHit_NHitsVsArm->Fill(1.,NHits[1]);
  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_PS_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_PS_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
