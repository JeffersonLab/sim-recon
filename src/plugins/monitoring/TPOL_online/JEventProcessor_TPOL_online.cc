#include <stdint.h>
#include <vector>

#include "JEventProcessor_TPOL_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseData.h>
#include <TRIGGER/DL1Trigger.h>
#include <TPOL/DTPOLSectorDigiHit.h>
#include <TPOL/DTPOLHit.h>
#include <TPOL/DTPOLHit_factory.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

// root hist pointers
const int Nsectors = DTPOLHit_factory::NSECTORS;
const int NmultBins = 25; // number of bins for multiplicity histograms
static TH1I *tpol_num_events;
static TH1I *hRaw_NHits;
static TH1I *hHit_NHits;
static TH1I *hHit_Occupancy;
static TH1F *hHit_Phi;
static TH1I *hHit_Peak;
static TH2I *hHit_PeakVsSector;
static TH1I *hHit_Integral;
static TH2I *hHit_IntegralVsSector;
static TH1I *hHit_Time;
static TH2I *hHit_TimeVsSector;
static TH2F *hHit_TimeVsPhi;
static TH2I *hHit_TimeVsPeak;
static TH2I *hHit_TimeVsIntegral;
//
static TH1I *hDigiHit_NHits;
static TH1I *hDigiHit_NSamplesPedestal;
static TH1I *hDigiHit_Pedestal;
static TProfile *hDigiHit_PedestalVsSector;
static TH1I *hDigiHit_QualityFactor;
static TH1I *hDigiHit_PulseNumber;
static TH2I *hDigiHit_PulseNumberVsSector;
static TH1I *hDigiHit_Occupancy;
static TH1I *hDigiHit_RawPeak;
static TH2I *hDigiHit_RawPeakVsSector;
static TH1I *hDigiHit_RawIntegral;
static TH2I *hDigiHit_RawIntegralVsSector;
static TH1I *hDigiHit_NSamplesIntegral;
static TH2I *hDigiHit_PeakVsSector;
static TH2I *hDigiHit_IntegralVsPeak;
static TH2I *hDigiHit_IntegralVsSector;
static TH1I *hDigiHit_PulseTime;
static TH1I *hDigiHit_Time;
static TH2I *hDigiHit_TimeVsSector;
static TH2I *hDigiHit_TimeVsPeak;
static TH2I *hDigiHit_TimeVsIntegral;
//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_TPOL_online());
    }
}


//----------------------------------------------------------------------------------


JEventProcessor_TPOL_online::JEventProcessor_TPOL_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TPOL_online::~JEventProcessor_TPOL_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TPOL_online::init(void) {

    // create root folder for TPOL and cd to it, store main dir
    TDirectory *mainDir = gDirectory;
    TDirectory *tpolDir = gDirectory->mkdir("TPOL");
    tpolDir->cd();
    // book hists
    tpol_num_events = new TH1I("tpol_num_events","TPOL number of events",1,0.5,1.5);
    hRaw_NHits = new TH1I("Raw_NHits",";TPOL raw hit multiplicity;raw hits;events",NmultBins,0.5,0.5+NmultBins);
    // hit-level hists (after calibration)
    gDirectory->mkdir("Hit")->cd();
    hHit_NHits = new TH1I("Hit_NHits","TPOL hit multiplicity;hits;events",NmultBins,0.5,0.5+NmultBins);
    hHit_Occupancy = new TH1I("Hit_Occupancy","TPOL occupancy;sector;hits / sector",Nsectors,0.5,0.5+Nsectors);
    hHit_Phi = new TH1F("Hit_Phi","TPOL azimuthal angle;triplet azimuthal angle [degrees];hits / sector",Nsectors,0.0,360.0);
    hHit_Peak = new TH1I("Hit_Peak","TPOL fADC pulse peak;pulse peak;hits",410,0.0,4100.0);
    hHit_PeakVsSector = new TH2I("Hit_PeakVsSector","TPOL fADC pulse peak vs. sector;sector;pulse peak",Nsectors,0.5,0.5+Nsectors,410,0.0,4100.0);
    hHit_Integral = new TH1I("Hit_Integral","TPOL fADC pulse integral;pulse integral;hits",1000,0.0,30000.0);
    hHit_IntegralVsSector = new TH2I("Hit_IntegralVsSector","TPOL fADC pulse integral vs. sector;sector;pulse integral",Nsectors,0.5,0.5+Nsectors,1000,0.0,30000.0);
    hHit_Time = new TH1I("Hit_Time","TPOL time;time [ns];hits / 2 ns",400,-400.0,400.0);
    hHit_TimeVsSector = new TH2I("Hit_TimeVsSector","TPOL time vs. sector;sector;time [ns]",Nsectors,0.5,0.5+Nsectors,400,-400.0,400.0);
    hHit_TimeVsPhi = new TH2F("Hit_TimeVsPhi","TPOL time vs. phi;#phi [degrees];time [ns]",Nsectors,0.0,360.0,400,-400.0,400.0);
    hHit_TimeVsIntegral = new TH2I("Hit_TimeVsIntegral","TPOL time vs. integral;pulse integral;time [ns]",500,0.0,30000.0,400,-400.0,400.0);
    hHit_TimeVsPeak = new TH2I("Hit_TimeVsPeak","TPOL time vs. peak;pulse peak;time [ns]",410,0.0,4100.0,400,-400.0,400.0);
    // digihit-level hists
    tpolDir->cd();
    gDirectory->mkdir("DigiHit")->cd();
    hDigiHit_NHits = new TH1I("DigiHit_NfadcHits","TPOL fADC hit multiplicity;raw hits;events",NmultBins,0.5,0.5+NmultBins);
    hDigiHit_NSamplesPedestal = new TH1I("DigiHit_NSamplesPedestal","TPOL fADC pedestal samples;pedestal samples;raw hits",50,-0.5,49.5);
    hDigiHit_Pedestal = new TH1I("DigiHit_Pedestal","TPOL fADC pedestals;pedestal [fADC counts];raw hits",200,0.0,200.0);
    hDigiHit_PedestalVsSector = new TProfile("DigiHit_PedestalVsSector","TPOL pedestal vs. sector;sector;average pedestal [fADC counts]",Nsectors,0.5,0.5+Nsectors,"s");
    hDigiHit_QualityFactor = new TH1I("DigiHit_QualityFactor","TPOL fADC quality factor;quality factor;raw hits",4,-0.5,3.5);
    hDigiHit_PulseNumber = new TH1I("DigiHit_PulseNumber","TPOL fADC pulse number;pulse number;raw hits",4,-0.5,3.5);
    hDigiHit_PulseNumberVsSector = new TH2I("DigiHit_PulseNumberVsSector","TPOL fADC pulse number vs. sector;sector;pulse number;raw hits",Nsectors,0.5,0.5+Nsectors,4,-0.5,3.5);
    hDigiHit_Occupancy = new TH1I("DigiHit_fadcOccupancy","TPOL fADC hit occupancy;sector;raw hits / counter",Nsectors,0.5,0.5+Nsectors);
    hDigiHit_RawPeak = new TH1I("DigiHit_RawPeak","TPOL fADC pulse peak (raw);pulse peak (raw);raw hits",410,0.0,4100.0);
    hDigiHit_RawPeakVsSector = new TH2I("DigiHit_RawPeakVsSector","TPOL fADC pulse peak (raw) vs. sector;sector;pulse peak (raw)",Nsectors,0.5,0.5+Nsectors,410,0.0,4100.0);
    hDigiHit_RawIntegral = new TH1I("DigiHit_RawIntegral","TPOL fADC pulse integral (raw);pulse integral (raw);raw hits",1000,0.0,30000.0);
    hDigiHit_RawIntegralVsSector = new TH2I("DigiHit_RawIntegralVsSector","TPOL fADC pulse integral (raw) vs. sector;sector;pulse integral (raw)",Nsectors,0.5,0.5+Nsectors,1000,0.0,30000.0);
    hDigiHit_NSamplesIntegral = new TH1I("DigiHit_NSamplesIntegral","TPOL fADC integral samples;integral samples;raw hits",60,-0.5,59.5);
    hDigiHit_PeakVsSector = new TH2I("DigiHit_PeakVsSector","TPOL fADC pulse peak vs. sector;sector;pulse peak",Nsectors,0.5,0.5+Nsectors,410,0.0,4100.0);
    hDigiHit_IntegralVsPeak = new TH2I("DigiHit_IntegralVsPeak","TPOL fADC pulse integral vs. peak;pulse peak;pulse integral",410,0.0,4100.0,1000,0.0,30000.0);
    hDigiHit_IntegralVsSector = new TH2I("DigiHit_IntegralVsSector","TPOL fADC pulse integral vs. sector;sector;pulse integral",Nsectors,0.5,0.5+Nsectors,1000,0.0,30000.0);
    hDigiHit_PulseTime = new TH1I("DigiHit_PulseTime","TPOL fADC pulse time;pulse time [62.5 ps];raw hits",1000,0.0,6500.0);
    hDigiHit_Time = new TH1I("DigiHit_Time","TPOL fADC pulse time;pulse time [ns];raw hits / 2 ns",200,0.0,400.0);
    hDigiHit_TimeVsSector = new TH2I("DigiHit_TimeVsSector","TPOL fADC pulse time vs. sector;sector;pulse time [ns]",Nsectors,0.5,0.5+Nsectors,200,0.0,400.0);
    hDigiHit_TimeVsPeak = new TH2I("DigiHit_TimeVsPeak","TPOL time vs. peak;pulse peak;time [ns]",410,0.0,4100.0,200,0.0,400.0);
    hDigiHit_TimeVsIntegral = new TH2I("DigiHit_fadcTimeVsIntegral","TPOL fADC pulse time vs. integral;pulse integral;pulse time [ns]",500,0.0,30000.0,200,0.0,400.0);
    // back to main dir
    mainDir->cd();

    return NOERROR;
}

//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
    // This is called whenever the run number changes
    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
    // This is called for every event. Use of common resources like writing
    // to a file or filling a histogram should be mutex protected. Using
    // loop-Get(...) to get reconstructed objects (and thereby activating the
    // reconstruction algorithm) should be done outside of any mutex lock
    // since multiple threads may call this method at the same time.
    const DL1Trigger *trig_words = NULL;
    uint32_t trig_mask, fp_trig_mask;
    try {
        eventLoop->GetSingle(trig_words);
    } catch(...) {};
    if (trig_words) {
        trig_mask = trig_words->trig_mask;
        fp_trig_mask = trig_words->fp_trig_mask;
    }
    else {
        trig_mask = 0;
        fp_trig_mask = 0;
    }
    int trig_bits = fp_trig_mask > 0 ? 10 + fp_trig_mask:trig_mask;
    // Select PS-triggered events
    if (trig_bits!=8) {
        return NOERROR;
    }
    vector<const DTPOLHit*> hits;
    eventLoop->Get(hits);
    vector<const DTPOLSectorDigiHit*> sdhits;
    eventLoop->Get(sdhits);
    vector<const Df250WindowRawData*> windowrawdata;
    eventLoop->Get(windowrawdata);
    // cache pulse pedestal and window raw data objects
    map< const DTPOLSectorDigiHit*, pair<const Df250PulsePedestal*, const Df250WindowRawData*> > pp_wr_cache;
    map< const DTPOLSectorDigiHit*, const Df250PulseData* > pd_cache;
    for(unsigned int i=0; i < sdhits.size(); i++) {
        const Df250PulsePedestal* pulsePed = NULL;
        const Df250PulseIntegral* pulseInt = NULL;
        const Df250PulseData* pulseDat = NULL;
        const Df250WindowRawData* rawData = NULL;
        sdhits[i]->GetSingle(pulsePed);
        sdhits[i]->GetSingle(pulseInt);
        sdhits[i]->GetSingle(pulseDat);
        if (pulseInt) pulseInt->GetSingle(rawData);
        pp_wr_cache[sdhits[i]] = pair<const Df250PulsePedestal*, const Df250WindowRawData*>(pulsePed, rawData);
        pd_cache[sdhits[i]] = pulseDat;
    }

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

    int NHits = 0;
    for(unsigned int i=0;i<windowrawdata.size();i++) {
        if(windowrawdata[i]->rocid == 84
            &&
            (windowrawdata[i]->slot == 13 || windowrawdata[i]->slot == 14)) {
                NHits++;
        }
    }
    hRaw_NHits->Fill(NHits);
    if(sdhits.size()>0) tpol_num_events->Fill(1);
    hDigiHit_NHits->Fill(sdhits.size());
    for(unsigned int i=0; i < sdhits.size(); i++) {
        //double ped = sdhits[i]->pedestal/sdhits[i]->nsamples_pedestal;
        hDigiHit_NSamplesPedestal->Fill(sdhits[i]->nsamples_pedestal);
        const Df250PulsePedestal* pulsePed = pp_wr_cache[sdhits[i]].first;
        const Df250PulseData* pulseDat = pd_cache[sdhits[i]];
        double ped = 0.0; double peak = 0.0;
        if (pulsePed) {
            ped = pulsePed->pedestal;
            peak = pulsePed->pulse_peak;
            hDigiHit_PulseNumber->Fill(pulsePed->pulse_number);
            hDigiHit_PulseNumberVsSector->Fill(sdhits[i]->sector,pulsePed->pulse_number);
        } else {
            ped = sdhits[i]->pedestal;
            peak = sdhits[i]->pulse_peak;
            if(pulseDat) {
                hDigiHit_PulseNumber->Fill(pulseDat->pulse_number);
                hDigiHit_PulseNumberVsSector->Fill(sdhits[i]->sector,pulseDat->pulse_number);
            }
        }
        hDigiHit_RawPeak->Fill(peak);
        hDigiHit_Pedestal->Fill(ped);
        if (ped==0.0||peak==0.0) continue;
        hDigiHit_PedestalVsSector->Fill(sdhits[i]->sector,ped);
        hDigiHit_Occupancy->Fill(sdhits[i]->sector);
        hDigiHit_RawPeakVsSector->Fill(sdhits[i]->sector,peak);
        hDigiHit_RawIntegral->Fill(sdhits[i]->pulse_integral);
        hDigiHit_RawIntegralVsSector->Fill(sdhits[i]->sector,sdhits[i]->pulse_integral);
        hDigiHit_NSamplesIntegral->Fill(sdhits[i]->nsamples_integral);
        double Pint = sdhits[i]->pulse_integral-sdhits[i]->nsamples_integral*ped; // pedestal-subtracted pulse integral
        hDigiHit_IntegralVsSector->Fill(sdhits[i]->sector,Pint);
        hDigiHit_IntegralVsPeak->Fill(peak-ped,Pint);
        hDigiHit_PeakVsSector->Fill(sdhits[i]->sector,peak-ped);
        hDigiHit_PulseTime->Fill(sdhits[i]->pulse_time);
        double t_ns = 0.0625*sdhits[i]->pulse_time;
        hDigiHit_Time->Fill(t_ns);
        hDigiHit_TimeVsSector->Fill(sdhits[i]->sector,t_ns);
        hDigiHit_TimeVsPeak->Fill(peak-ped,t_ns);
        hDigiHit_TimeVsIntegral->Fill(Pint,t_ns);
        hDigiHit_QualityFactor->Fill(sdhits[i]->QF);
    }
    hHit_NHits->Fill(hits.size());
    for(unsigned int i=0; i<hits.size(); i++) {
        hHit_Occupancy->Fill(hits[i]->sector);
        hHit_Phi->Fill(hits[i]->phi);
        hHit_Integral->Fill(hits[i]->integral);
        hHit_IntegralVsSector->Fill(hits[i]->sector,hits[i]->integral);
        hHit_Peak->Fill(hits[i]->pulse_peak);
        hHit_PeakVsSector->Fill(hits[i]->sector,hits[i]->pulse_peak);
        hHit_Time->Fill(hits[i]->t);
        hHit_TimeVsSector->Fill(hits[i]->sector,hits[i]->t);
        hHit_TimeVsPhi->Fill(hits[i]->phi,hits[i]->t);
        hHit_TimeVsIntegral->Fill(hits[i]->integral,hits[i]->t);
        hHit_TimeVsPeak->Fill(hits[i]->pulse_peak,hits[i]->t);
    }

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}
//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::erun(void) {
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::fini(void) {
    // Called before program exit after event processing is finished.
    return NOERROR;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
