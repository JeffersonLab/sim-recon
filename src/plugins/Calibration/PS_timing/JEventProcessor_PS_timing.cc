// $Id$
//
//    File: JEventProcessor_PS_timing.cc
// Created: Sat Nov 21 17:21:28 EST 2015
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.el7.x86_64 x86_64)
//

#include "JEventProcessor_PS_timing.h"

using namespace std;
using namespace jana;

#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <PAIR_SPECTROMETER/DPSPair.h>
#include <TAGGER/DTAGHHit.h>
#include <RF/DRFTime.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>

const int NC_PSC = 16; // number of PSC modules (counters)
const int NC_PS = 290; // number of PS columns (tiles)
const int NC_TAGH = 274; // number of TAGH counter slots

static TH2I *hPSC_tdcadcTimeDiffVsID;
static TH2I *hPSCRF_tdcTimeDiffVsID;
static TH2I *hPSRF_adcTimeDiffVsID;
static TH2I *hTAGHRF_tdcTimeDiffVsID;

extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_PS_timing());
    }
} // "C"


//------------------
// JEventProcessor_PS_timing (Constructor)
//------------------
JEventProcessor_PS_timing::JEventProcessor_PS_timing()
{

}

//------------------
// ~JEventProcessor_PS_timing (Destructor)
//------------------
JEventProcessor_PS_timing::~JEventProcessor_PS_timing()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_PS_timing::init(void)
{
    // This is called once at program startup. If you are creating
    // and filling historgrams in this plugin, you should lock the
    // ROOT mutex like this:

    const double Tl = -200.0;
    const double Th = 600.0;
    const int NTb = 8000;
    TDirectory *mainDir = gDirectory;
    TDirectory *psDir = gDirectory->mkdir("PS_timing");
    psDir->cd();
    hPSC_tdcadcTimeDiffVsID = new TH2I("PSC_tdcadcTimeDiffVsID","PSC TDC-ADC time difference vs. counter ID;counter ID;time(TDC) - time(ADC) [ns]",NC_PSC,0.5,0.5+NC_PSC,NTb,Tl,Th);
    hPSCRF_tdcTimeDiffVsID = new TH2I("PSCRF_tdcTimeDiffVsID","PSC-RF TDC time difference vs. counter ID;counter ID;time(TDC) - time(RF) [ns]",NC_PSC,0.5,0.5+NC_PSC,NTb,Tl,Th);
    hPSRF_adcTimeDiffVsID = new TH2I("PSRF_adcTimeDiffVsID","PS-RF ADC time difference vs. counter ID;counter ID;time(ADC) - time(RF) [ns]",NC_PS,0.5,0.5+NC_PS,NTb,Tl,Th);
    hTAGHRF_tdcTimeDiffVsID = new TH2I("TAGHRF_tdcTimeDiffVsID","TAGH-RF TDC time difference vs. counter ID;counter ID;time(TDC) - time(RF) [ns]",NC_TAGH,0.5,0.5+NC_TAGH,NTb,Tl,Th);
    mainDir->cd();

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_PS_timing::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_PS_timing::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // This is called for every event. Use of common resources like writing
    // to a file or filling a histogram should be mutex protected. Using
    // loop->Get(...) to get reconstructed objects (and thereby activating the
    // reconstruction algorithm) should be done outside of any mutex lock
    // since multiple threads may call this method at the same time.
    vector<const DPSCPair*> cpairs;
    loop->Get(cpairs);
    vector<const DPSPair*> fpairs;
    loop->Get(fpairs);

    vector<const DTAGHHit*> taghhits;
    loop->Get(taghhits, "Calib");

    const DRFTime* rfTime = nullptr;
    vector <const DRFTime*> rfTimes;
    loop->Get(rfTimes, "PSC");
    if (rfTimes.size() > 0)
        rfTime = rfTimes[0];
    else
        return NOERROR;

    // FILL HISTOGRAMS
    // Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
    japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

    double t_RF = rfTime->dTime;
    if (cpairs.size() >= 1) { // PSC
        const DPSCHit* clhit = cpairs[0]->ee.first; // left hit in coarse PS
        const DPSCHit* crhit = cpairs[0]->ee.second;// right hit in coarse PS
        hPSC_tdcadcTimeDiffVsID->Fill(clhit->module,clhit->t-clhit->time_fadc);
        hPSCRF_tdcTimeDiffVsID->Fill(clhit->module,clhit->t-t_RF);
        hPSC_tdcadcTimeDiffVsID->Fill(crhit->module+8,crhit->t-crhit->time_fadc);
        hPSCRF_tdcTimeDiffVsID->Fill(crhit->module+8,crhit->t-t_RF);
        for (const auto& h : taghhits) hTAGHRF_tdcTimeDiffVsID->Fill(h->counter_id,h->t-t_RF);
        if (fpairs.size() >= 1) { // PS
            const DPSHit* flhit = fpairs[0]->ee.first;  // left hit in fine PS
            const DPSHit* frhit = fpairs[0]->ee.second; // right hit in fine PS
            hPSRF_adcTimeDiffVsID->Fill(flhit->column,flhit->t-t_RF);
            hPSRF_adcTimeDiffVsID->Fill(frhit->column+145,frhit->t-t_RF);
        }
    }

    japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_PS_timing::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_PS_timing::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

