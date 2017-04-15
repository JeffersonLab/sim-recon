// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>

#include "JEventProcessor_TAGH_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "TTAB/DTTabUtilities.h"
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGHGeometry.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/DEPICSvalue.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

// Nslots: total number of TAGH counter slots
// 274 total (218 have counters for default GlueX configuration)
const int Nslots = DTAGHGeometry::kCounterCount;
// counter (slot) id =  HV id for hv id from 1-131 (upstream of microscope)
// counter (slot) id =  41 + HV id for hv id from 132-233 (downstream of microscope)
const int NHVchannels = 233;
const int NupstreamCounters = 131;
//
const int NmultBins = 300; // number of bins for multiplicity histograms
//
double beam_current = -1.0;
const double counts_cut = 1000.0;
// root hist pointers
static TH1I *tagh_num_events;
static TH1I *hBeamCurrent;
static TH2F *hHit_HasTDCvsHasADC;
static TH1I *hHit_RawNHits;
static TH1I *hHit_NHits;
static TH1I *hHit_NHits_us;
static TH1I *hHit_NHits_ds;
static TH1F *hHit_Occupancy;
static TH2I *hHit_HVidVsSlotID;
static TH1F *hHit_Energy;
static TH2I *hHit_EnergyVsSlotID;
static TH1I *hHit_Integral;
static TH2I *hHit_IntegralVsSlotID;
//static TH1I *hHit_fadcNpe;
static TH1I *hHit_fadcTime;
static TH2I *hHit_fadcTimeVsSlotID;
static TH1I *hHit_Time;
static TH2F *hHit_TimeVsSlotID;
static TH2F *hHit_TimeVsEnergy;
static TH2I *hHit_TimeVsIntegral;
static TH1I *hHit_tdcTime;
static TH2I *hHit_tdcTimeVsSlotID;
static TH2I *hHit_tdcadcTimeDiffVsSlotID;
static TH2I *hHit_tdcadcTimeDiffVsIntegral;
//
static TH1I *hDigiHit_NfadcHits;
static TH2I *hDigiHit_NfadcHitsVsSlotID;
static TH1I *hDigiHit_NfadcHits_multiPeak;
static TH1I *hDigiHit_NSamplesPedestal;
static TH1I *hDigiHit_Pedestal;
static TProfile *hDigiHit_PedestalVsSlotID;
static TH1I *hDigiHit_QualityFactor;
static TH1I *hDigiHit_PulseNumber;
static TH2I *hDigiHit_PulseNumberVsSlotID;
static TH1I *hDigiHit_fadcOccupancy;
static TH1I *hDigiHit_RawPeak;
static TH2I *hDigiHit_RawPeakVsSlotID;
static TH1I *hDigiHit_RawIntegral;
static TH2I *hDigiHit_RawIntegralVsSlotID;
static TH1I *hDigiHit_NSamplesIntegral;
static TH2I *hDigiHit_PeakVsSlotID;
static TH2I *hDigiHit_IntegralVsPeak;
static TH2I *hDigiHit_IntegralVsSlotID;
static TH1I *hDigiHit_PulseTime;
static TH1I *hDigiHit_fadcTime;
static TH2I *hDigiHit_fadcTimeVsSlotID;
static TH2I *hDigiHit_fadcTimeVsIntegral;
static TH1I *hDigiHit_NtdcHits;
static TH2I *hDigiHit_NtdcHitsVsSlotID;
static TH1I *hDigiHit_tdcOccupancy;
static TH1I *hDigiHit_tdcRawTime;
static TH1I *hDigiHit_tdcTime;
static TH2I *hDigiHit_tdcTimeVsSlotID;
static TH2I *hDigiHit_tdcTimeVsfadcTime;
static TH1I *hDigiHit_tdcadcTimeDiff;
static TH2I *hDigiHit_tdcadcTimeDiffVsSlotID;
static TH2I *hDigiHit_tdcadcTimeDiffVsIntegral;
//
static TH1I *hDigiHit_NfadcHits_cut;
static TH1I *hDigiHit_fadcOccupancy_cut;
static TH1I *hDigiHit_fadcTime_cut;
static TH2I *hDigiHit_fadcTimeVsSlotID_cut;
static TH2I *hDigiHit_fadcTimeVsQF_cut;
static TH2I *hDigiHit_adctdcMatchesVsSlotID_cut;
//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_TAGH_online());
    }
}


//----------------------------------------------------------------------------------


JEventProcessor_TAGH_online::JEventProcessor_TAGH_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TAGH_online::~JEventProcessor_TAGH_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TAGH_online::init(void) {

    // create root folder for tagh and cd to it, store main dir
    TDirectory *mainDir = gDirectory;
    TDirectory *taghDir = gDirectory->mkdir("TAGH");
    taghDir->cd();

    // book hists
    tagh_num_events = new TH1I("tagh_num_events","TAGH number of events",1,0.5,1.5);
    hBeamCurrent = new TH1I("BeamCurrent","Beam current;beam current [nA]",600,0.0,300.0);
    // hit-level hists (after calibration)
    gDirectory->mkdir("Hit")->cd();
    hHit_HasTDCvsHasADC = new TH2F("Hit_HasTDCvsHasADC","TAGH has TDC? vs. has ADC?;fADC status;TDC status",2,-0.5,1.5,2,-0.5,1.5);
    hHit_RawNHits = new TH1I("Hit_RawNHits","TAGH hit multiplicity (raw);hits (raw);events",NmultBins,0.5,0.5+NmultBins);
    //
    hHit_NHits = new TH1I("Hit_NHits","TAGH hit multiplicity;hits;events",NmultBins,0.5,0.5+NmultBins);
    hHit_NHits_us = new TH1I("Hit_NHits_us","TAGH upstream hit multiplicity;upstream hits;events",NmultBins,0.5,0.5+NmultBins);
    hHit_NHits_ds = new TH1I("Hit_NHits_ds","TAGH downstream hit multiplicity;downstream hits;events",NmultBins,0.5,0.5+NmultBins);
    hHit_Occupancy    = new TH1F("Hit_Occupancy","TAGH occupancy;counter (slot) ID;hits / counter",Nslots,0.5,0.5+Nslots);
    hHit_HVidVsSlotID    = new TH2I("Hit_HVidVsSlotID","TAGH high voltage channel ID vs. counter ID;counter (slot) ID;high voltage channel ID",Nslots,0.5,0.5+Nslots,NHVchannels,0.5,0.5+NHVchannels);
    hHit_Energy    = new TH1F("Hit_Energy","TAGH energy;photon energy [GeV];hits / counter",120,0.0,12.0);
    hHit_EnergyVsSlotID    = new TH2I("Hit_EnergyVsSlotID","TAGH energy vs. counter ID;counter (slot) ID;photon energy [GeV]",Nslots,0.5,0.5+Nslots,120,0.0,12.0);
    hHit_Integral    = new TH1I("Hit_Integral","TAGH fADC pulse integral;pulse integral;hits",1000,0.0,30000.0);
    hHit_IntegralVsSlotID    = new TH2I("Hit_IntegralVsSlotID","TAGH fADC pulse integral vs. counter ID;counter (slot) ID;pulse integral",Nslots,0.5,0.5+Nslots,1000,0.0,30000.0);
    //hHit_fadcNpe    = new TH1I("Hit_fadcNpe","TAGH fADC number of photoelectrons;photoelectrons;hits",1000,0.0,30000.0);
    hHit_fadcTime    = new TH1I("Hit_fadcTime","TAGH fADC time;time [ns];hits / 400 ps",2000,-400.0,400.0);
    hHit_fadcTimeVsSlotID    = new TH2I("Hit_fadcTimeVsSlotID","TAGH fADC time vs. counter ID;counter (slot) ID;time [ns]",Nslots,0.5,0.5+Nslots,2000,-400.0,400.0);
    hHit_Time    = new TH1I("Hit_Time","TAGH time;time [ns];hits / 400 ps",2000,-400.0,400.0);
    hHit_TimeVsSlotID    = new TH2F("Hit_TimeVsSlotID","TAGH time vs. counter ID;counter (slot) ID;time [ns]",Nslots,0.5,0.5+Nslots,2000,-400.0,400.0);
    hHit_TimeVsEnergy    = new TH2F("Hit_TimeVsEnergy","TAGH time vs. energy;energy [GeV];time [ns]",120,0.0,12.0,2000,-400.0,400.0);
    hHit_TimeVsIntegral    = new TH2I("Hit_TimeVsIntegral","TAGH time vs. integral;pulse integral;time [ns]",500,0.0,15000.0,2000,-400.0,400.0);
    hHit_tdcTime    = new TH1I("Hit_tdcTime","TAGH TDC time;time [ns];hits / 400 ps",2000,-400.0,400.0);
    hHit_tdcTimeVsSlotID    = new TH2I("Hit_tdcTimeVsSlotID","TAGH TDC time vs. counter ID;counter (slot) ID;time [ns]",Nslots,0.5,0.5+Nslots,2000,-400.0,400.0);
    hHit_tdcadcTimeDiffVsSlotID    = new TH2I("Hit_tdcadcTimeDiffVsSlotID","TAGH TDC/ADC time difference vs. counter ID;counter (slot) ID;time(TDC) - time(ADC) [ns]",Nslots,0.5,0.5+Nslots,200,-40.0,40.0);
    hHit_tdcadcTimeDiffVsIntegral    = new TH2I("Hit_tdcadcTimeDiffVsIntegral","TAGH TDC/ADC time difference vs. integral;pulse integral;time(TDC) - time(ADC) [ns]",500,0.0,15000.0,200,-40.0,40.0);
    // digihit-level hists
    taghDir->cd();
    gDirectory->mkdir("DigiHit")->cd();
    hDigiHit_NfadcHits = new TH1I("DigiHit_NfadcHits","TAGH fADC hit multiplicity;raw hits;events",NmultBins,0.5,0.5+NmultBins);
    hDigiHit_NSamplesPedestal    = new TH1I("DigiHit_NSamplesPedestal","TAGH fADC pedestal samples;pedestal samples;raw hits",50,-0.5,49.5);
    hDigiHit_Pedestal = new TH1I("DigiHit_Pedestal","TAGH fADC pedestals;pedestal [fADC counts];raw hits",200,0.0,200.0);
    hDigiHit_PedestalVsSlotID = new TProfile("DigiHit_PedestalVsSlotID","TAGH pedestal vs. counter ID;counter (slot) ID;average pedestal [fADC counts]",Nslots,0.5,0.5+Nslots,"s");
    hDigiHit_QualityFactor = new TH1I("DigiHit_QualityFactor","TAGH fADC quality factor;quality factor;raw hits",4,-0.5,3.5);
    hDigiHit_PulseNumber = new TH1I("DigiHit_PulseNumber","TAGH fADC pulse number;pulse number;raw hits",4,-0.5,3.5);
    hDigiHit_PulseNumberVsSlotID = new TH2I("DigiHit_PulseNumberVsSlotID","TAGH fADC pulse number vs. counter ID;counter (slot) ID;pulse number;raw hits",Nslots,0.5,0.5+Nslots,4,-0.5,3.5);
    hDigiHit_fadcOccupancy    = new TH1I("DigiHit_fadcOccupancy","TAGH fADC hit occupancy;counter (slot) ID;raw hits / counter",Nslots,0.5,0.5+Nslots);
    hDigiHit_RawPeak    = new TH1I("DigiHit_RawPeak","TAGH fADC pulse peak (raw);pulse peak (raw);raw hits",410,0.0,4100.0);
    hDigiHit_RawPeakVsSlotID    = new TH2I("DigiHit_RawPeakVsSlotID","TAGH fADC pulse peak (raw) vs. counter ID;counter (slot) ID;pulse peak (raw)",Nslots,0.5,0.5+Nslots,410,0.0,4100.0);
    hDigiHit_RawIntegral    = new TH1I("DigiHit_RawIntegral","TAGH fADC pulse integral (raw);pulse integral (raw);raw hits",1000,0.0,30000.0);
    hDigiHit_RawIntegralVsSlotID    = new TH2I("DigiHit_RawIntegralVsSlotID","TAGH fADC pulse integral (raw) vs. counter ID;counter (slot) ID;pulse integral (raw)",Nslots,0.5,0.5+Nslots,1000,0.0,30000.0);
    hDigiHit_NSamplesIntegral    = new TH1I("DigiHit_NSamplesIntegral","TAGH fADC integral samples;integral samples;raw hits",60,-0.5,59.5);
    hDigiHit_PeakVsSlotID    = new TH2I("DigiHit_PeakVsSlotID","TAGH fADC pulse peak vs. counter ID;counter (slot) ID;pulse peak",Nslots,0.5,0.5+Nslots,410,0.0,4100.0);
    hDigiHit_IntegralVsPeak    = new TH2I("DigiHit_IntegralVsPeak","TAGH fADC pulse integral vs. peak;pulse peak;pulse integral",410,0.0,4100.0,1000,0.0,30000.0);
    hDigiHit_IntegralVsSlotID    = new TH2I("DigiHit_IntegralVsSlotID","TAGH fADC pulse integral vs. counter ID;counter (slot) ID;pulse integral",Nslots,0.5,0.5+Nslots,1000,0.0,30000.0);
    hDigiHit_PulseTime    = new TH1I("DigiHit_PulseTime","TAGH fADC pulse time;pulse time [62.5 ps];raw hits",1000,0.0,6500.0);
    hDigiHit_fadcTime    = new TH1I("DigiHit_fadcTime","TAGH fADC pulse time;pulse time [ns];raw hits / 400 ps",1000,0.0,400.0);
    hDigiHit_fadcTimeVsSlotID    = new TH2I("DigiHit_fadcTimeVsSlotID","TAGH fADC pulse time vs. counter ID;counter (slot) ID;pulse time [ns]",Nslots,0.5,0.5+Nslots,1000,0.0,400.0);
    hDigiHit_fadcTimeVsIntegral    = new TH2I("DigiHit_fadcTimeVsIntegral","TAGH fADC pulse time vs. integral;pulse integral;pulse time [ns]",500,0.0,15000.0,1000,0.0,400.0);
    hDigiHit_NtdcHits = new TH1I("DigiHit_NtdcHits","TAGH TDC hit multiplicity;raw hits;events",NmultBins,0.5,0.5+NmultBins);
    hDigiHit_NtdcHitsVsSlotID = new TH2I("DigiHit_NtdcHitsVsSlotID","TAGH TDC hit multiplicity vs. counter ID;counter ID;raw hits",Nslots,0.5,0.5+Nslots,8,0.5,8.5);
    hDigiHit_tdcOccupancy    = new TH1I("DigiHit_tdcOccupancy","TAGH TDC hit occupancy;counter (slot) ID;raw hits / counter",Nslots,0.5,0.5+Nslots);
    hDigiHit_tdcRawTime    = new TH1I("DigiHit_tdcRawTime","TAGH TDC raw time;time [60 ps];raw hits",1000,0.0,65500.0);
    hDigiHit_tdcTime    = new TH1I("DigiHit_tdcTime","TAGH TDC time;time [ns];raw hits / 400 ps",2000,0.0,800.0);
    hDigiHit_tdcTimeVsSlotID    = new TH2I("DigiHit_tdcTimeVsSlotID","TAGH TDC time vs. counter ID;counter ID;TDC time [ns]",Nslots,0.5,0.5+Nslots,2000,0.0,800.0);
    hDigiHit_tdcTimeVsfadcTime    = new TH2I("DigiHit_tdcTimeVsfadcTime","TAGH TDC time vs. ADC time;fADC time [ns];TDC time [ns]",400,0.0,400.0,800,0.0,800.0);
    hDigiHit_tdcadcTimeDiff    = new TH1I("DigiHit_tdcadcTimeDiff","TAGH TDC/ADC time difference;time(TDC) - time(ADC) [ns];raw hits / 400 ps",1000,0.0,400.0);
    hDigiHit_tdcadcTimeDiffVsSlotID    = new TH2I("DigiHit_tdcadcTimeDiffVsSlotID","TAGH TDC/ADC time difference vs. counter ID;counter ID;time(TDC) - time(ADC) [ns]",Nslots,0.5,0.5+Nslots,1000,0.0,400.0);
    hDigiHit_tdcadcTimeDiffVsIntegral    = new TH2I("DigiHit_tdcadcTimeDiffVsIntegral","TAGH TDC/ADC time difference vs. pulse integral;pulse integral;time(TDC) - time(ADC) [ns]",500,0.0,15000.0,1000,0.0,400.0);
    // allow multiple peaks in window and cut on ADC counts
    hDigiHit_NfadcHits_multiPeak = new TH1I("DigiHit_NfadcHits_multiPeak","TAGH fADC hit multiplicity (> 400 counts,multiple peaks allowed);raw hits;events",NmultBins,0.5,0.5+NmultBins);
    hDigiHit_NfadcHitsVsSlotID = new TH2I("DigiHit_NfadcHitsVsSlotID","TAGH fADC hit multiplicity vs. counter ID (> 400 counts);counter ID;raw hits",Nslots,0.5,0.5+Nslots,8,0.5,8.5);
    // digihit-level hists after cut on ADC integral counts
    hDigiHit_NfadcHits_cut = new TH1I("DigiHit_NfadcHits_cut","TAGH fADC hit multiplicity (> 1k ADC integral counts);raw hits;events",NmultBins,0.5,0.5+NmultBins);
    hDigiHit_fadcOccupancy_cut    = new TH1I("DigiHit_fadcOccupancy_cut","TAGH fADC hit occupancy (> 1k ADC integral counts);counter (slot) ID;raw hits / counter",Nslots,0.5,0.5+Nslots);
    hDigiHit_fadcTime_cut    = new TH1I("DigiHit_fadcTime_cut","TAGH fADC pulse time (> 1k ADC integral counts);pulse time [ns];raw hits / 400 ps",1000,0.0,400.0);
    hDigiHit_fadcTimeVsSlotID_cut    = new TH2I("DigiHit_fadcTimeVsSlotID_cut","TAGH fADC pulse time vs. counter ID (> 1k ADC integral counts);counter (slot) ID;pulse time [ns]",Nslots,0.5,0.5+Nslots,1000,0.0,400.0);
    hDigiHit_fadcTimeVsQF_cut    = new TH2I("DigiHit_fadcTimeVsQF_cut","TAGH fADC pulse time vs. quality factor (> 1k ADC integral counts);fADC quality factor;pulse time [ns]",4,-0.5,3.5,1000,0.0,400.0);
    hDigiHit_adctdcMatchesVsSlotID_cut = new TH2I("DigiHit_adctdcMatchesVsSlotID_cut","TAGH #TDC matches / fADC hit vs. counter ID (> 1k ADC integral counts);counter ID;#TDC matches / fADC hit",Nslots,0.5,0.5+Nslots,8,-0.5,7.5);
    // back to main dir
    mainDir->cd();

    return NOERROR;
}

//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGH_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
    // This is called whenever the run number changes
    // extract the TAGH geometry
    vector<const DTAGHGeometry*> taghGeomVect;
    eventLoop->Get(taghGeomVect);
    if (taghGeomVect.size() == 0) return OBJECT_NOT_AVAILABLE;
    const DTAGHGeometry& taghGeom = *(taghGeomVect[0]);
    // get photon energy bin low of each counter for energy histogram binning
    double Elows[Nslots + 1];
    double slots[Nslots + 1];
    for (int i = 0; i < Nslots; i++) {
        Elows[i] = taghGeom.getElow(Nslots - i);
        slots[i] = 0.5 + i;
    }
    // add the upper limit
    Elows[Nslots] = taghGeom.getEhigh(1);
    slots[Nslots] = 0.5 + Nslots;

    const int Ntime = 2000;
    double Tlows[Ntime + 1];
    for (int i = 0; i <= Ntime; i++) {
        Tlows[i] = -400.0 + i*0.4;
    }

    // FILL HISTOGRAMS
    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

    if (hHit_Energy->GetEntries() == 0) hHit_Energy->SetBins(Nslots,Elows);
    if (hHit_EnergyVsSlotID->GetEntries() == 0) hHit_EnergyVsSlotID->SetBins(Nslots,slots,Nslots,Elows);
    if (hHit_TimeVsEnergy->GetEntries() == 0) hHit_TimeVsEnergy->SetBins(Nslots,Elows,Ntime,Tlows);

    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGH_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
    // This is called for every event. Use of common resources like writing
    // to a file or filling a histogram should be mutex protected. Using
    // loop-Get(...) to get reconstructed objects (and thereby activating the
    // reconstruction algorithm) should be done outside of any mutex lock
    // since multiple threads may call this method at the same time.

    // Get TAGH hits and digihits
    vector<const DTAGHHit*> hits;
    eventLoop->Get(hits, "Calib");
    vector<const DTAGHDigiHit*> digihits;
    eventLoop->Get(digihits);
    vector<const DTAGHTDCDigiHit*> tdcdigihits;
    eventLoop->Get(tdcdigihits);

    const DTTabUtilities* ttabUtilities = nullptr;
    eventLoop->GetSingle(ttabUtilities);

    // Cache pulse data and window raw data objects
    map< const DTAGHDigiHit*, pair<const Df250PulseData*, const Df250WindowRawData*> > pd_wrd_cache;
    for (const auto& hit : digihits) {
        const Df250PulseData* pd = nullptr;
        const Df250WindowRawData* wrd = nullptr;
        hit->GetSingle(pd);
        if (pd != nullptr) pd->GetSingle(wrd);
        pd_wrd_cache[hit] = pair<const Df250PulseData*, const Df250WindowRawData*>(pd, wrd);
    }

    // Get beam current
    vector<const DEPICSvalue*> epicsVals;
    eventLoop->Get(epicsVals);
    for (const auto& ev : epicsVals) {
        if (ev->name == "IBCAD00CRCUR6") {
            beam_current = ev->fval;
            break;
        }
    }

    // FILL HISTOGRAMS
    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

    hBeamCurrent->Fill(beam_current);
    if (digihits.size() > 0 || tdcdigihits.size() > 0)
        tagh_num_events->Fill(1);

    hHit_RawNHits->Fill(hits.size());
    hDigiHit_NfadcHits->Fill(digihits.size());
    hDigiHit_NtdcHits->Fill(tdcdigihits.size());
    int NfadcHits_cut = 0;
    int NHits_hasADC = 0;  int NHits_hasADC_us = 0;  int NHits_hasADC_ds = 0;
    int Nadc[Nslots];
    for (int i = 0; i < Nslots; i++) Nadc[i] = 0;
    for (const auto& hit : digihits) {
        double ped = (double)hit->pedestal/hit->nsamples_pedestal;
        hDigiHit_NSamplesPedestal->Fill(hit->nsamples_pedestal);
        hDigiHit_Pedestal->Fill(ped);
        hDigiHit_RawPeak->Fill(hit->pulse_peak);
        if (ped == 0.0 || hit->pulse_peak == 0) continue;
        hDigiHit_PedestalVsSlotID->Fill(hit->counter_id,ped);
        hDigiHit_fadcOccupancy->Fill(hit->counter_id);
        hDigiHit_RawPeakVsSlotID->Fill(hit->counter_id,hit->pulse_peak);
        hDigiHit_RawIntegral->Fill(hit->pulse_integral);
        hDigiHit_RawIntegralVsSlotID->Fill(hit->counter_id,hit->pulse_integral);
        hDigiHit_NSamplesIntegral->Fill(hit->nsamples_integral);
        double pI = hit->pulse_integral-hit->nsamples_integral*ped;
        hDigiHit_IntegralVsSlotID->Fill(hit->counter_id,pI);
        hDigiHit_IntegralVsPeak->Fill(hit->pulse_peak-ped,pI);
        hDigiHit_PeakVsSlotID->Fill(hit->counter_id,hit->pulse_peak-ped);
        hDigiHit_PulseTime->Fill(hit->pulse_time);
        double t_ns = 0.0625*hit->pulse_time;
        hDigiHit_fadcTime->Fill(t_ns);
        hDigiHit_fadcTimeVsSlotID->Fill(hit->counter_id,t_ns);
        hDigiHit_fadcTimeVsIntegral->Fill(pI,t_ns);
        hDigiHit_QualityFactor->Fill(hit->QF);
        const Df250PulseData* pd = pd_wrd_cache[hit].first;
        int pN = (pd != nullptr) ? pd->pulse_number : -1;
        if (pN > -1) hDigiHit_PulseNumber->Fill(pN);
        if (pN > -1) hDigiHit_PulseNumberVsSlotID->Fill(hit->counter_id,pN);
        if (pI > counts_cut)  {
            NfadcHits_cut++;
            hDigiHit_fadcOccupancy_cut->Fill(hit->counter_id);
            hDigiHit_fadcTime_cut->Fill(t_ns);
            hDigiHit_fadcTimeVsSlotID_cut->Fill(hit->counter_id,t_ns);
            hDigiHit_fadcTimeVsQF_cut->Fill(hit->QF,t_ns);
        }
        const Df250WindowRawData* wrd = pd_wrd_cache[hit].second;
        size_t threshold = 400;
        if (wrd != nullptr) {
            size_t prevSample = 100;
            int Ncrosses = 0;
            if (!wrd->invalid_samples && !wrd->overflow) {
                for (size_t sample : wrd->samples) {
                    if (sample > threshold && prevSample < threshold) {
                        Ncrosses++;
                    }
                    prevSample = sample;
                }
            }
            Nadc[hit->counter_id-1] = Ncrosses;
        }
        else {
            if (hit->pulse_peak > threshold) Nadc[hit->counter_id-1]++;
        }
    }

    int NmultiPeak = 0;
    if (digihits.size() > 0) {
        for (int i = 0; i < Nslots; i++) {
            NmultiPeak += Nadc[i];
            hDigiHit_NfadcHitsVsSlotID->Fill(i+1,Nadc[i]);
        }
    }
    hDigiHit_NfadcHits_multiPeak->Fill(NmultiPeak);
    hDigiHit_NfadcHits_cut->Fill(NfadcHits_cut);
    for (const auto& hit : digihits) {
        double ped = (double)hit->pedestal/hit->nsamples_pedestal;
        double pI = hit->pulse_integral-ped*hit->nsamples_integral;
        if (hit->pedestal > 0 && pI > counts_cut) {
            int matches = 0;
            for (const auto& tdchit : tdcdigihits) {
                if (hit->counter_id == tdchit->counter_id) {
                    matches++;
                    double T_tdc = ttabUtilities->Convert_DigiTimeToNs_F1TDC(tdchit);
                    double T_adc = 0.0625*hit->pulse_time;
                    hDigiHit_tdcTimeVsfadcTime->Fill(T_adc,T_tdc);
                    hDigiHit_tdcadcTimeDiff->Fill(T_tdc-T_adc);
                    hDigiHit_tdcadcTimeDiffVsSlotID->Fill(hit->counter_id,T_tdc-T_adc);
                    hDigiHit_tdcadcTimeDiffVsIntegral->Fill(pI,T_tdc-T_adc);
                }
            }
            hDigiHit_adctdcMatchesVsSlotID_cut->Fill(hit->counter_id,matches);
        }
    }

    int Ntdc[Nslots];
    for (int i = 0; i < Nslots; i++) Ntdc[i] = 0;
    for (const auto& hit : tdcdigihits) {
        Ntdc[hit->counter_id-1]++;
        hDigiHit_tdcOccupancy->Fill(hit->counter_id);
        hDigiHit_tdcRawTime->Fill(hit->time);
        double T_tdc = ttabUtilities->Convert_DigiTimeToNs_F1TDC(hit);
        hDigiHit_tdcTime->Fill(T_tdc);
        hDigiHit_tdcTimeVsSlotID->Fill(hit->counter_id,T_tdc);
    }
    if (tdcdigihits.size() > 0) {
        for (int i = 0; i < Nslots; i++) {
            hDigiHit_NtdcHitsVsSlotID->Fill(i+1,Ntdc[i]);
        }
    }

    // Fill histograms with calibrated-hit data
    for (const auto& hit : hits) {
        hHit_HasTDCvsHasADC->Fill(hit->has_fADC,hit->has_TDC);
        if (hit->has_fADC) {
            NHits_hasADC++;
            if (hit->counter_id<=NupstreamCounters) NHits_hasADC_us++;
            else NHits_hasADC_ds++;
            hHit_Occupancy->Fill(hit->counter_id);
            int HVid = (hit->counter_id <= NupstreamCounters) ? hit->counter_id : (hit->counter_id - Nslots + NHVchannels);
            hHit_HVidVsSlotID->Fill(hit->counter_id,HVid);
            hHit_Energy->Fill(hit->E);
            hHit_EnergyVsSlotID->Fill(hit->counter_id,hit->E);
            hHit_Integral->Fill(hit->integral);
            hHit_IntegralVsSlotID->Fill(hit->counter_id,hit->integral);
            //hHit_fadcNpe->Fill(hit->npe_fadc);
            hHit_fadcTime->Fill(hit->time_fadc);
            hHit_fadcTimeVsSlotID->Fill(hit->counter_id,hit->time_fadc);
            hHit_Time->Fill(hit->t);
            hHit_TimeVsSlotID->Fill(hit->counter_id,hit->t);
            hHit_TimeVsEnergy->Fill(hit->E,hit->t);
            hHit_TimeVsIntegral->Fill(hit->integral,hit->t);
            if (hit->has_TDC) {
                hHit_tdcTime->Fill(hit->time_tdc);
                hHit_tdcTimeVsSlotID->Fill(hit->counter_id,hit->time_tdc);
                hHit_tdcadcTimeDiffVsSlotID->Fill(hit->counter_id,hit->time_tdc-hit->time_fadc);
                hHit_tdcadcTimeDiffVsIntegral->Fill(hit->integral,hit->time_tdc-hit->time_fadc);
            }
        }
    }
    hHit_NHits->Fill(NHits_hasADC);
    hHit_NHits_us->Fill(NHits_hasADC_us);
    hHit_NHits_ds->Fill(NHits_hasADC_ds);

    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGH_online::erun(void) {
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGH_online::fini(void) {
    // Called before program exit after event processing is finished.
    return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
