// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>
#include <TMath.h>


#include "JEventProcessor_TAGM_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>

// Define some constants
//const uint32_t NROWS = 5;
const uint32_t NCOLUMNS = 102;
const uint32_t NSINGLES = 20;

const float MIN_ADC_PINT_LOG10 = 0.;
const float MAX_ADC_PINT_LOG10 = 5.;
const uint32_t BINCOUNT_ADC_PINT = 200;
//const float ADC_PINT_PER_PIXEL = 6.5;
const float MIN_HIT_NPIX = 0.;
const float MAX_HIT_NPIX = 1000.;
const uint32_t BINCOUNT_HIT_NPIX = 200;

const float MIN_ADC_TIME = 0.;
const float MAX_ADC_TIME = 16000.;
const uint32_t BINCOUNT_ADC_TIME = 16000;
const float_t ADC_NS_PER_COUNT = 0.0625;
const float MIN_HIT_TADC = 0.;
const float MAX_HIT_TADC = 1000.;
const uint32_t BINCOUNT_HIT_TADC = 16000;
const float MIN_HIT_TIME = 0.;
const float MAX_HIT_TIME = 4000.;
const uint32_t BINCOUNT_HIT_TIME = 65000;

const float MIN_TDC_TIME = 0.;
const float MAX_TDC_TIME = 65000.;
const uint32_t BINCOUNT_TDC_TIME = 65000;
const float_t TDC_NS_PER_COUNT = 0.060;

const float MIN_ADC_PED = 0.;
const float MAX_ADC_PED = 800.;
const uint32_t BINCOUNT_ADC_PED = 800;

const float MIN_ADC_QF = 0.;
const float MAX_ADC_QF = 4.;
const uint32_t BINCOUNT_ADC_QF = 5;

const float MIN_ADC_NSI = 0.;
const float MAX_ADC_NSI = 150.;
const uint32_t BINCOUNT_ADC_NSI = 151;

const float MIN_ADC_NSP = 0.;
const float MAX_ADC_NSP = 20.;
const uint32_t BINCOUNT_ADC_NSP = 21;

// DTAGDigiHit histogram pointers
static TH1F *tagm_adc_seen,     *tagms_adc_seen;
static TH1F *tagm_adc_pint,     *tagms_adc_pint;
static TH2F *tagm_adc_pint_2d,  *tagms_adc_pint_2d;
static TH1F *tagm_adc_time,     *tagms_adc_time;
static TH2F *tagm_adc_time_2d,  *tagms_adc_time_2d;
static TH1F *tagm_adc_ped,      *tagms_adc_ped;
static TH2F *tagm_adc_ped_2d,   *tagms_adc_ped_2d;
static TH1F *tagm_adc_qf,       *tagms_adc_qf;
static TH2F *tagm_adc_qf_2d,    *tagms_adc_qf_2d;
static TH1F *tagm_adc_nsi,      *tagms_adc_nsi;
static TH2F *tagm_adc_nsi_2d,   *tagms_adc_nsi_2d;
static TH1F *tagm_adc_nsp,      *tagms_adc_nsp;
static TH2F *tagm_adc_nsp_2d,   *tagms_adc_nsp_2d;
static TH2F *tagm_adc_pint_nsi, *tagms_adc_pint_nsi;
static TH2F *tagm_adc_ped_nsp,  *tagms_adc_ped_nsp;

// DSCTDCDigiHit histogram pointers
static TH1F *tagm_tdc_seen,     *tagms_tdc_seen;
static TH1F *tagm_tdc_time,     *tagms_tdc_time;
static TH2F *tagm_tdc_time_2d,  *tagms_tdc_time_2d;

// DTAGMHit histogram pointers
static TH1F *tagm_hit_seen,     *tagms_hit_seen;
static TH1F *tagm_hit_npix,     *tagms_hit_npix;
static TH2F *tagm_hit_npix_2d,  *tagms_hit_npix_2d;
static TH1F *tagm_hit_tadc,     *tagms_hit_tadc;
static TH2F *tagm_hit_tadc_2d,  *tagms_hit_tadc_2d;
static TH1F *tagm_hit_time,     *tagms_hit_time;
static TH2F *tagm_hit_time_2d,  *tagms_hit_time_2d;
static TH2F *tagm_hit_time_tadc, *tagms_hit_time_tadc;

// Dynamic arrays of histogram pointers
TH1F** tagm_adc_pint_col  = new TH1F*[NCOLUMNS];
TH1F** tagm_adc_time_col  = new TH1F*[NCOLUMNS];
TH1F** tagm_adc_ped_col   = new TH1F*[NCOLUMNS];
TH1F** tagm_adc_qf_col    = new TH1F*[NCOLUMNS];
TH1F** tagm_tdc_time_col  = new TH1F*[NCOLUMNS];
TH1F** tagm_adc_mult_col  = new TH1F*[NCOLUMNS];
TH1F** tagm_tdc_mult_col  = new TH1F*[NCOLUMNS];
TH1F** tagm_hit_time_col  = new TH1F*[NCOLUMNS];

TH1F** tagms_adc_pint_sng  = new TH1F*[NSINGLES];
TH1F** tagms_adc_time_sng  = new TH1F*[NSINGLES];
TH1F** tagms_adc_ped_sng   = new TH1F*[NSINGLES];
TH1F** tagms_adc_qf_sng    = new TH1F*[NSINGLES];
TH1F** tagms_tdc_time_sng  = new TH1F*[NSINGLES];
TH1F** tagms_adc_mult_sng  = new TH1F*[NSINGLES];
TH1F** tagms_tdc_mult_sng  = new TH1F*[NSINGLES];
TH1F** tagms_hit_time_sng  = new TH1F*[NSINGLES];

// Number of events histogram pointer
static TH1F *tagm_num_events;

// Time difference histograms
static TH1F *tagm_hit_tdiff,     *tagms_hit_tdiff;
static TH1F *tagm_tdc_tdiff_all, *tagms_tdc_tdiff_all;

// Multiplicity histograms
static TH1F *tagm_adc_mult;
static TH1F *tagm_tdc_mult;
static TH2F *tagm_adc_mult_2d,  *tagms_adc_mult_2d;
static TH2F *tagm_tdc_mult_2d,  *tagms_tdc_mult_2d;
static TH2F *tagm_adc_tdc_mult;


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TAGM_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_TAGM_online::JEventProcessor_TAGM_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TAGM_online::~JEventProcessor_TAGM_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TAGM_online::init(void) {

  // lock all root operations
  japp->RootWriteLock();


  // create root folder for tagm and cd to it, store main dir
  TDirectory *main = gDirectory;
  TDirectory *tagmdir = gDirectory->mkdir("tagm");
  tagmdir->cd();
 
  // Book DTAGMDigiHit histosgrams 
  tagm_adc_seen    = new TH1F("tagm_adc_seen", 
                              "TAGM FADC250 column occupancy",
                              NCOLUMNS, 0., NCOLUMNS + 1.);
  tagm_adc_pint    = new TH1F("tagm_adc_pint",
                              "TAGM FADC250 pulse integral (log10)",
                              BINCOUNT_ADC_PINT, 
                              MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagm_adc_pint_2d = new TH2F("tagm_adc_pint_2d",
                              "TAGM FADC250 pulse integral (log10) vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_ADC_PINT, 
                              MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagm_adc_time    = new TH1F("tagm_adc_time",
                              "TAGM FADC250 pulse time counter",
                              BINCOUNT_ADC_TIME, MIN_ADC_TIME, MAX_ADC_TIME);
  tagm_adc_time_2d = new TH2F("tagm_adc_time_2d",
                              "TAGM FADC250 pulse time counter vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_ADC_TIME, MIN_ADC_TIME, MAX_ADC_TIME);
  tagm_adc_ped     = new TH1F("tagm_adc_ped",
                              "TAGM FADC250 pulse pedestal",
                              BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);
  tagm_adc_ped_2d  = new TH2F("tagm_adc_ped_2d",
                              "TAGM FADC250 pulse pedestal vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);
  tagm_adc_qf      = new TH1F("tagm_adc_qf", 
                              "TAGM FADC250 quality factor",
                              BINCOUNT_ADC_QF, 
                              MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
  tagm_adc_qf_2d   = new TH2F("tagm_adc_qf_2d", 
                              "TAGM FADC250 quality factor vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_ADC_QF, 
                              MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
  tagm_adc_nsi     = new TH1F("tagm_adc_nsi",
                              "TAGM number of samples used in pulse integral",
                              BINCOUNT_ADC_NSI, 
                              MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5);
  tagm_adc_nsi_2d  = new TH2F("tagm_adc_nsi_2d",
                              "TAGM samples used in pulse integral vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.0,
                              BINCOUNT_ADC_NSI, 
                              MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5);
  tagm_adc_nsp     = new TH1F("tagm_adc_nsp",
                              "TAGM number of samples used in pulse pedestal",
                              BINCOUNT_ADC_NSP, 
                              MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5);
  tagm_adc_nsp_2d  = new TH2F("tagm_adc_nsp_2d",
                              "TAGM samples used in pulse pedestal vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_ADC_NSP, 
                              MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5);
  tagm_adc_pint_nsi = new TH2F("tagm_adc_pint_nsi",
                              "TAGM pulse integral (log10) vs number of samples",
                              BINCOUNT_ADC_NSI, 
                              MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5,
                              BINCOUNT_ADC_PINT, 
                              MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagm_adc_ped_nsp = new TH2F("tagm_adc_ped_nsp",
                              "TAGM pulse pedestal vs pedestal samples",
                              BINCOUNT_ADC_NSP, 
                              MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5,
                              BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);

  // Book DTAGMTDCDigiHit histosgrams 
  tagm_tdc_seen    = new TH1F("tagm_tdc_seen", 
                              "TAGM F1TDC column occupancy",
                              NCOLUMNS, 0., NCOLUMNS + 1.);
  tagm_tdc_time    = new TH1F("tagm_tdc_time",
                              "TAGM F1TDC pulse time counter",
                              BINCOUNT_TDC_TIME, MIN_TDC_TIME, MAX_TDC_TIME);
  tagm_tdc_time_2d = new TH2F("tagm_tdc_time_2d",
                              "TAGM F1TDC pulse time counter vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_TDC_TIME, MIN_TDC_TIME, MAX_TDC_TIME);

  // Book DTAGMHit histosgrams
  tagm_hit_seen    = new TH1F("tagm_hit_seen", 
                              "TAGM hit column occupancy", 
                              NCOLUMNS, 0., NCOLUMNS + 1.);
  tagm_hit_npix    = new TH1F("tagm_hit_npix",
                              "TAGM hit pulse height (pixels)",
                              BINCOUNT_HIT_NPIX, MIN_HIT_NPIX, MAX_HIT_NPIX);
  tagm_hit_npix_2d = new TH2F("tagm_hit_npix_2d",
                              "TAGM hit pulse height (pixels) vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_HIT_NPIX, MIN_HIT_NPIX, MAX_HIT_NPIX);
  tagm_hit_tadc    = new TH1F("tagm_hit_tadc",
                              "TAGM hit pulse time (ns)",
                              BINCOUNT_HIT_TADC, MIN_HIT_TADC, MAX_HIT_TADC);
  tagm_hit_tadc_2d = new TH2F("tagm_hit_tadc_2d",
                              "TAGM hit pulse time (ns) vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_HIT_TADC, MIN_HIT_TADC, MAX_HIT_TADC);
  tagm_hit_time    = new TH1F("tagm_hit_time",
                              "TAGM hit discriminator time (ns)",
                              BINCOUNT_HIT_TIME, MIN_HIT_TIME, MAX_HIT_TIME);
  tagm_hit_time_2d = new TH2F("tagm_hit_time_2d",
                              "TAGM hit discriminator time (ns) vs column",
                              NCOLUMNS, 0., NCOLUMNS + 1.,
                              BINCOUNT_HIT_TIME, MIN_HIT_TIME, MAX_HIT_TIME);
  tagm_hit_time_tadc = new TH2F("tagm_hit_time_tadc",
                              "TAGM pulse time vs discriminator time (ns)",
                              1000, MIN_HIT_TIME, MAX_HIT_TIME,
                              1000, MIN_HIT_TADC, MAX_HIT_TADC);

  // Book dynamic arrays of histograms
  TDirectory *adc_pint_dir = gDirectory->mkdir("tagm_adc_pint_col");
  TDirectory *adc_time_dir = gDirectory->mkdir("tagm_adc_time_col");
  TDirectory *adc_ped_dir  = gDirectory->mkdir("tagm_adc_ped_col");
  TDirectory *adc_qf_dir   = gDirectory->mkdir("tagm_adc_qf_col");
  TDirectory *tdc_time_dir = gDirectory->mkdir("tagm_tdc_time_col");
  TDirectory *adc_mult_dir = gDirectory->mkdir("tagm_adc_mult_col");
  TDirectory *tdc_mult_dir = gDirectory->mkdir("tagm_tdc_mult_col");
  TDirectory *hit_time_dir = gDirectory->mkdir("tagm_hit_time_col");
  for (unsigned int i = 0; i < NCOLUMNS; i++) {
    adc_pint_dir->cd();
    tagm_adc_pint_col[i] = new TH1F(Form("tagm_adc_pint_col_%i", i+1),
                                    Form("FADC250 pulse integral (log10)"
                                         " for TAGM column %i", i+1),
                                    BINCOUNT_ADC_PINT, 
                                    MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
    adc_time_dir->cd();
    tagm_adc_time_col[i] = new TH1F(Form("tagm_adc_time_col_%i", i+1),
                                    Form("FADC250 pulse time counter"
                                         " for TAGM column %i", i+1),
                                    BINCOUNT_ADC_TIME, 
                                    MIN_ADC_TIME, MAX_ADC_TIME);
    adc_ped_dir->cd();
    tagm_adc_ped_col[i]  = new TH1F(Form("tagm_adc_ped_col_%i", i+1),
                                    Form("FADC250 pulse pedestal"
                                         " for TAGM column %i", i+1),
                                    BINCOUNT_ADC_PED, 
                                    MIN_ADC_PED, MAX_ADC_PED);
    adc_qf_dir->cd();
    tagm_adc_qf_col[i]   = new TH1F(Form("tagm_adc_qf_col_%i", i+1),
                                    Form("FADC250 pulse QF"
                                         " for TAGM column %i", i+1),
                                    BINCOUNT_ADC_QF, 
                                    MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
    tdc_time_dir->cd();
    tagm_tdc_time_col[i] = new TH1F(Form("tagm_tdc_time_col_%i", i+1),
                                    Form("F1TDC hit time counter"
                                         " for TAGM column %i", i+1),
                                    BINCOUNT_TDC_TIME, 
                                    MIN_TDC_TIME, MAX_TDC_TIME);
    adc_mult_dir->cd();
    tagm_adc_mult_col[i] = new TH1F(Form("tagm_adc_mult_col_%i", i+1),
                                    Form("FADC250 pulse multiplicity"
                                         " for TAGM column %i", i+1),
                                    10, 0., 10.);
    tdc_mult_dir->cd();
    tagm_tdc_mult_col[i] = new TH1F(Form("tagm_tdc_mult_col_%i", i+1),
                                    Form("F1TDC hit multiplicity"
                                         " for TAGM column %i", i+1),
                                    10, 0., 10.);
    hit_time_dir->cd();
    tagm_hit_time_col[i] = new TH1F(Form("tagm_hit_time_col_%i", i+1),
                                    Form("hit time for TAGM column %i", i+1),
                                    BINCOUNT_ADC_TIME, 
                                    MIN_ADC_TIME, MAX_ADC_TIME);
  }
  tagmdir->cd();

  // Book time difference histograms
  tagm_hit_tdiff     = new TH1F("tagm_hit_tdiff",
                                "TAGM TDC-ADC hit time difference, nearby pairs",
                                2*BINCOUNT_HIT_TIME/100,
                                -MAX_HIT_TIME/100., MAX_HIT_TIME/100.);
  tagm_tdc_tdiff_all = new TH1F("tagm_tdc_tdiff_all",
                                "TAGM F1TDC time difference, all pairs",
                                2*BINCOUNT_HIT_TIME,
                                -MAX_HIT_TIME, MAX_HIT_TIME);

  // Book multiplicity histograms
  tagm_tdc_mult      = new TH1F("tagm_tdc_mult",
                                "TAGM F1TDC multiplicity",
                                10, 0.5, 10.5);
  tagm_adc_mult      = new TH1F("tagm_adc_mult",
                                "TAGM FADC250 multiplicity",
                                10, 0.5, 10.5);
  tagm_adc_mult_2d   = new TH2F("tagm_adc_mult_2d",
                                "TAGM FADC250 multiplicity vs column",
                                NCOLUMNS, 0., NCOLUMNS + 1., 10, 0.5, 10.5);
  tagm_tdc_mult_2d   = new TH2F("tagm_tdc_mult_2d",
                                "TAGM F1TDC multiplicity vs column",
                                NCOLUMNS, 0., NCOLUMNS + 1., 10, 0.5, 10.5);
  tagm_adc_tdc_mult  = new TH2F("tagm_adc_tdc_mult",
                                "TAGM FADC250 vs. F1TDC multiplicity",
                                10, 0.5, 10.5, 10, 0.5, 10.5);

  // Book total number of events histogram
  tagm_num_events = new TH1F("tagm_num_events",
                            "TAGM number of events", 1, 0.5, 1.5);

  // Book histograms for individual fiber channels
  TDirectory *singles_dir = gDirectory->mkdir("singles (low-gain only)");
  singles_dir->cd();
  tagms_adc_seen    = new TH1F("tagms_adc_seen", 
                               "TAGM FADC250 singles occupancy",
                               NSINGLES, 0., NSINGLES + 1.);
  tagms_adc_pint    = new TH1F("tagms_adc_pint",
                               "TAGM FADC250 pulse integral (log10)",
                               BINCOUNT_ADC_PINT, 
                               MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagms_adc_pint_2d = new TH2F("tagms_adc_pint_2d",
                               "TAGM FADC250 pulse integral (log10) vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_ADC_PINT, 
                               MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagms_adc_time    = new TH1F("tagms_adc_time",
                               "TAGM FADC250 pulse time counter",
                               BINCOUNT_ADC_TIME, MIN_ADC_TIME, MAX_ADC_TIME);
  tagms_adc_time_2d = new TH2F("tagms_adc_time_2d",
                               "TAGM FADC250 pulse time counter vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_ADC_TIME, MIN_ADC_TIME, MAX_ADC_TIME);
  tagms_adc_ped     = new TH1F("tagms_adc_ped",
                               "TAGM FADC250 pulse pedestal",
                               BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);
  tagms_adc_ped_2d  = new TH2F("tagms_adc_ped_2d",
                               "TAGM FADC250 pulse pedestal vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);
  tagms_adc_qf      = new TH1F("tagms_adc_qf", 
                               "TAGM FADC250 quality factor",
                               BINCOUNT_ADC_QF, 
                               MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
  tagms_adc_qf_2d   = new TH2F("tagms_adc_qf_2d", 
                               "TAGM FADC250 quality factor vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_ADC_QF, 
                               MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
  tagms_adc_nsi     = new TH1F("tagms_adc_nsi",
                               "TAGM number of samples used in pulse integral",
                               BINCOUNT_ADC_NSI, 
                               MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5);
  tagms_adc_nsi_2d  = new TH2F("tagms_adc_nsi_2d",
                               "TAGM samples used in pulse integral vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.0,
                               BINCOUNT_ADC_NSI, 
                               MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5);
  tagms_adc_nsp     = new TH1F("tagms_adc_nsp",
                               "TAGM number of samples used in pulse pedestal",
                               BINCOUNT_ADC_NSP, 
                               MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5);
  tagms_adc_nsp_2d  = new TH2F("tagms_adc_nsp_2d",
                               "TAGM samples used in pulse pedestal vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_ADC_NSP, 
                               MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5);
  tagms_adc_pint_nsi = new TH2F("tagms_adc_pint_nsi",
                               "TAGM pulse integral (log10) vs number of samples",
                               BINCOUNT_ADC_NSI, 
                               MIN_ADC_NSI - 0.5, MAX_ADC_NSI + 0.5,
                               BINCOUNT_ADC_PINT, 
                               MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
  tagms_adc_ped_nsp = new TH2F("tagms_adc_ped_nsp",
                               "TAGM pulse pedestal vs pedestal samples",
                               BINCOUNT_ADC_NSP, 
                               MIN_ADC_NSP - 0.5, MAX_ADC_NSP + 0.5,
                               BINCOUNT_ADC_PED, MIN_ADC_PED, MAX_ADC_PED);

  // Book DTAGMTDCDigiHit histosgrams 
  tagms_tdc_seen    = new TH1F("tagms_tdc_seen", 
                               "TAGM F1TDC column occupancy",
                               NSINGLES, 0., NSINGLES + 1.);
  tagms_tdc_time    = new TH1F("tagms_tdc_time",
                               "TAGM F1TDC pulse time counter",
                               BINCOUNT_TDC_TIME, MIN_TDC_TIME, MAX_TDC_TIME);
  tagms_tdc_time_2d = new TH2F("tagms_tdc_time_2d",
                               "TAGM F1TDC pulse time counter vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_TDC_TIME, MIN_TDC_TIME, MAX_TDC_TIME);

  // Book DTAGMHit histosgrams
  tagms_hit_seen    = new TH1F("tagms_hit_seen", 
                               "TAGM hit fiber no. occupancy", 
                               NSINGLES, 0., NSINGLES + 1.);
  tagms_hit_npix    = new TH1F("tagms_hit_npix",
                               "TAGM hit pulse height (pixels)",
                               BINCOUNT_HIT_NPIX, MIN_HIT_NPIX, MAX_HIT_NPIX);
  tagms_hit_npix_2d = new TH2F("tagms_hit_npix_2d",
                               "TAGM hit pulse height (pixels) vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_HIT_NPIX, MIN_HIT_NPIX, MAX_HIT_NPIX);
  tagms_hit_tadc    = new TH1F("tagms_hit_tadc",
                               "TAGM hit pulse time (ns)",
                               BINCOUNT_HIT_TADC, MIN_HIT_TADC, MAX_HIT_TADC);
  tagms_hit_tadc_2d = new TH2F("tagms_hit_tadc_2d",
                               "TAGM hit pulse time (ns) vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_HIT_TADC, MIN_HIT_TADC, MAX_HIT_TADC);
  tagms_hit_time    = new TH1F("tagms_hit_time",
                               "TAGM hit discriminator time (ns)",
                               BINCOUNT_HIT_TIME, MIN_HIT_TIME, MAX_HIT_TIME);
  tagms_hit_time_2d = new TH2F("tagms_hit_time_2d",
                               "TAGM hit discriminator time (ns) vs fiber no.",
                               NSINGLES, 0., NSINGLES + 1.,
                               BINCOUNT_HIT_TIME, MIN_HIT_TIME, MAX_HIT_TIME);
  tagms_hit_time_tadc = new TH2F("tagms_hit_time_tadc",
                               "TAGM pulse time vs discriminator time (ns)",
                               1000, MIN_HIT_TIME, MAX_HIT_TIME,
                               1000, MIN_HIT_TADC, MAX_HIT_TADC);

  // Book dynamic arrays of histograms
  TDirectory *adc_pint_sdir = gDirectory->mkdir("tagms_adc_pint_sng");
  TDirectory *adc_time_sdir = gDirectory->mkdir("tagms_adc_time_sng");
  TDirectory *adc_ped_sdir  = gDirectory->mkdir("tagms_adc_ped_sng");
  TDirectory *adc_qf_sdir   = gDirectory->mkdir("tagms_adc_qf_sng");
  TDirectory *tdc_time_sdir = gDirectory->mkdir("tagms_tdc_time_sng");
  TDirectory *adc_mult_sdir = gDirectory->mkdir("tagms_adc_mult_sng");
  TDirectory *tdc_mult_sdir = gDirectory->mkdir("tagms_tdc_mult_sng");
  TDirectory *hit_time_sdir = gDirectory->mkdir("tagms_hit_time_sng");
  for (unsigned int i = 0; i < NSINGLES; i++) {
    adc_pint_sdir->cd();
    tagms_adc_pint_sng[i] = new TH1F(Form("tagms_adc_pint_sng_%i", i+1),
                                     Form("FADC250 pulse integral (log10)"
                                          " for TAGM fiber no. %i", i+1),
                                     BINCOUNT_ADC_PINT, 
                                     MIN_ADC_PINT_LOG10, MAX_ADC_PINT_LOG10);
    adc_time_sdir->cd();
    tagms_adc_time_sng[i] = new TH1F(Form("tagms_adc_time_sng_%i", i+1),
                                     Form("FADC250 pulse time counter"
                                          " for TAGM fiber no. %i", i+1),
                                     BINCOUNT_ADC_TIME, 
                                     MIN_ADC_TIME, MAX_ADC_TIME);
    adc_ped_sdir->cd();
    tagms_adc_ped_sng[i]  = new TH1F(Form("tagms_adc_ped_sng_%i", i+1),
                                     Form("FADC250 pulse pedestal"
                                          " for TAGM fiber no. %i", i+1),
                                     BINCOUNT_ADC_PED, 
                                     MIN_ADC_PED, MAX_ADC_PED);
    adc_qf_sdir->cd();
    tagms_adc_qf_sng[i]   = new TH1F(Form("tagms_adc_qf_sng_%i", i+1),
                                     Form("FADC250 pulse QF"
                                          " for TAGM fiber no. %i", i+1),
                                     BINCOUNT_ADC_QF, 
                                     MIN_ADC_QF - 0.5, MAX_ADC_QF + 0.5);
    tdc_time_sdir->cd();
    tagms_tdc_time_sng[i] = new TH1F(Form("tagms_tdc_time_sng_%i", i+1),
                                     Form("F1TDC hit time counter"
                                          " for TAGM fiber no. %i", i+1),
                                     BINCOUNT_TDC_TIME, 
                                     MIN_TDC_TIME, MAX_TDC_TIME);
    adc_mult_sdir->cd();
    tagms_adc_mult_sng[i] = new TH1F(Form("tagms_adc_mult_sng_%i", i+1),
                                     Form("FADC250 pulse multiplicity"
                                          " for TAGM fiber no. %i", i+1),
                                     10, 0., 10.);
    tdc_mult_sdir->cd();
    tagms_tdc_mult_sng[i] = new TH1F(Form("tagms_tdc_mult_sng_%i", i+1),
                                     Form("F1TDC hit multiplicity"
                                          " for TAGM fiber no. %i", i+1),
                                     10, 0., 10.);
    hit_time_sdir->cd();
    tagms_hit_time_sng[i] = new TH1F(Form("tagms_hit_time_sng_%i", i+1),
                                     Form("hit time for TAGM fiber no. %i", i+1),
                                     BINCOUNT_ADC_TIME, 
                                     MIN_ADC_TIME, MAX_ADC_TIME);
  }
  singles_dir->cd();

  // Book time difference histograms
  tagms_hit_tdiff     = new TH1F("tagms_hit_tdiff",
                                 "TAGM TDC-ADC hit time difference, nearby pairs",
                                 2*BINCOUNT_HIT_TIME/100,
                                 -MAX_HIT_TIME/100., MAX_HIT_TIME/100.);
  tagms_tdc_tdiff_all = new TH1F("tagms_tdc_tdiff_all",
                                 "TAGM F1TDC time difference, all pairs",
                                 2*BINCOUNT_HIT_TIME,
                                 -MAX_HIT_TIME, MAX_HIT_TIME);

  // Book multiplicity histograms
  tagms_adc_mult_2d   = new TH2F("tagms_adc_mult_2d",
                                 "TAGM FADC250 multiplicity vs fiber no.",
                                 NSINGLES, 0., NSINGLES + 1., 10, 0.5, 10.5);
  tagms_tdc_mult_2d   = new TH2F("tagms_tdc_mult_2d",
                                 "TAGM F1TDC multiplicity vs fiber no.",
                              NSINGLES, 0., NSINGLES + 1., 10, 0.5, 10.5);

  // back to main dir
  main->cd();

  // unlock
  japp->RootUnLock();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGM_online::brun(JEventLoop *eventLoop, int runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGM_online::evnt(JEventLoop *eventLoop, int eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  std::vector<const DTAGMDigiHit*> digihits;             // FADC250 DigiHits
  std::vector<const DTAGMTDCDigiHit*> tdcdigihits;       // F1TDC DigiHits
  std::vector<const DTAGMHit*> hits;                     // hits
  eventLoop->Get(digihits);
  eventLoop->Get(tdcdigihits);
  eventLoop->Get(hits);

  // Lock ROOT mutex so other threads won't interfere 
  japp->RootWriteLock();

  // Histogram the total number of events
  if( (digihits.size()>0) || (tdcdigihits.size()>0) )
	  tagm_num_events->Fill(1);

  // Histogram the ADC and TDC multiplicities 
  tagm_adc_mult->Fill(digihits.size());
  tagm_tdc_mult->Fill(tdcdigihits.size());
  tagm_adc_tdc_mult->Fill(tdcdigihits.size(), digihits.size());

  // Prepare multiplicity counters
  int column_adc_hits[NCOLUMNS];
  int column_tdc_hits[NCOLUMNS];
  for (unsigned int i=0; i < NCOLUMNS; ++i)
    column_adc_hits[i] = column_tdc_hits[i] = 0;
  int column_adc_hits_total = 0;
  int column_tdc_hits_total = 0;

  // DTAGMDigiHits
  std::vector<const DTAGMDigiHit*>::iterator iter;
  for (iter = digihits.begin(); iter != digihits.end(); ++iter) {
    int row = (*iter)->row;
    int column = (*iter)->column;
    uint32_t pedestal_avg = (*iter)->pedestal / (*iter)->nsamples_pedestal;
    if (row == 0) {

      // Fill 1D Histograms
      tagm_adc_seen->Fill((*iter)->column);
      tagm_adc_pint->Fill(log10((*iter)->pulse_integral));
      tagm_adc_time->Fill((*iter)->pulse_time);
      tagm_adc_ped->Fill(pedestal_avg);
      tagm_adc_qf->Fill((*iter)->QF);
      tagm_adc_nsi->Fill((*iter)->nsamples_integral);
      tagm_adc_nsp->Fill((*iter)->nsamples_pedestal);

      // Fill 2D histograms
      tagm_adc_pint_2d->Fill(column, log10((*iter)->pulse_integral));
      tagm_adc_time_2d->Fill(column, (*iter)->pulse_time);
      tagm_adc_ped_2d->Fill(column, pedestal_avg);
      tagm_adc_qf_2d->Fill(column, (*iter)->QF);
      tagm_adc_nsi_2d->Fill(column, (*iter)->nsamples_integral);
      tagm_adc_nsp_2d->Fill(column, (*iter)->nsamples_pedestal);
      tagm_adc_pint_nsi->Fill((*iter)->nsamples_integral,
                              log10((*iter)->pulse_integral));
      tagm_adc_ped_nsp->Fill((*iter)->nsamples_pedestal,
                             pedestal_avg);

      // Fill dynamic array of histograms
      tagm_adc_pint_col[column - 1]->Fill(log10((*iter)->pulse_integral));
      tagm_adc_time_col[column - 1]->Fill((*iter)->pulse_time);
      tagm_adc_ped_col[column - 1]->Fill(pedestal_avg);
      tagm_adc_qf_col[column - 1]->Fill((*iter)->QF);

      // Calculate the FADC250 multiplicities for each channel
      ++column_adc_hits[column - 1];
      ++column_adc_hits_total;
    }
  }

  // DSCTDCDigiHits
  std::vector<const DTAGMTDCDigiHit*>::iterator titer;
  for (titer = tdcdigihits.begin(); titer != tdcdigihits.end(); ++titer) {
    int row = (*titer)->row;
    int column = (*titer)->column;
    if (row == 0) {

      // Fill 1D histograms
      tagm_tdc_seen->Fill((*titer)->column);
      tagm_tdc_time->Fill((*titer)->time);

      // Fill 2D histograms
      tagm_tdc_time_2d->Fill((*titer)->column, (*titer)->time);

      // Fill dynamic array of 1D histograms
      tagm_tdc_time_col[column - 1]->Fill((*titer)->time);

      // Calculate the F1TDC multiplicities for each channel
      ++column_tdc_hits[column - 1];
      ++column_tdc_hits_total;
    }
  }

  // Fill the 2D multiplicity histograms for ADC and TDC
  for (unsigned int col = 0; col < NCOLUMNS; ++col){
    if (column_adc_hits[col] > 0) {
      tagm_adc_mult_2d->Fill(col + 1, column_adc_hits[col]);
      tagm_adc_mult_col[col]->Fill(column_adc_hits[col]);
    }

    if (column_tdc_hits[col] > 0) {
      tagm_tdc_mult_2d->Fill(col + 1, column_tdc_hits[col]);
      tagm_tdc_mult_col[col]->Fill(column_tdc_hits[col]);
    }
  }

  // DTAGMHits
  std::vector<const DTAGMHit*>::iterator hiter;
  for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
    int row = (*hiter)->row;
    int column = (*hiter)->column;
    if (row == 0) {

      // Fill 1D histograms
      tagm_hit_seen->Fill((*hiter)->column);
      tagm_hit_npix->Fill(log10((*hiter)->npix_fadc));
      tagm_hit_tadc->Fill((*hiter)->time_fadc);
      tagm_hit_time->Fill((*hiter)->t);
      tagm_hit_time_tadc->Fill((*hiter)->t, (*hiter)->time_fadc);
      tagm_hit_tdiff->Fill((*hiter)->t - (*hiter)->time_fadc);
 
      // Fill 2D histograms
      tagm_hit_npix_2d->Fill((*hiter)->column, (*hiter)->npix_fadc);
      tagm_hit_tadc_2d->Fill((*hiter)->column, (*hiter)->time_fadc);
      tagm_hit_time_2d->Fill((*hiter)->column, (*hiter)->t);

      // Fill dynamic array of 1D histograms
      tagm_hit_time_col[column - 1]->Fill((*hiter)->t);
    }
  }

  // Find all possible time differences for multiple TDC hits
  if ((column_tdc_hits_total > 1) && (column_adc_hits_total > 1) ) {
    for (titer = tdcdigihits.begin(); titer != tdcdigihits.end(); ++titer) {
      for (iter = digihits.begin(); iter != digihits.end(); ++iter) {
        if ((*iter)->column == (*titer)->column &&
            (*iter)->row == 0 && (*titer)->row ==0) 
        {
          tagm_tdc_tdiff_all->Fill((*titer)->time * TDC_NS_PER_COUNT -
                                   (*iter)->pulse_time * ADC_NS_PER_COUNT);
        }
      }
    }
  }

  // Now repeat all of the above for single fiber channels

  // Prepare multiplicity counters
  int single_adc_hits[NSINGLES];
  int single_tdc_hits[NSINGLES];
  for (unsigned int i=0; i < NSINGLES; ++i)
    single_adc_hits[i] = single_tdc_hits[i] = 0;
  int single_adc_hits_total = 0;
  int single_tdc_hits_total = 0;

  // DTAGMDigiHits
  for (iter = digihits.begin(); iter != digihits.end(); ++iter) {
    int row = (*iter)->row;
    int column = (*iter)->column;
    uint32_t pedestal_avg = (*iter)->pedestal / (*iter)->nsamples_pedestal;
    int fiberNo = (row < 1)?      0 :
                  (column == 7)?  row :
                  (column == 25)? row + 5 :
                  (column == 79)? row + 10 :
                  (column == 97)? row + 15 : 0;
    if (fiberNo > 0) {

      // Fill 1D Histograms
      tagms_adc_seen->Fill(fiberNo);
      tagms_adc_pint->Fill(log10((*iter)->pulse_integral));
      tagms_adc_time->Fill((*iter)->pulse_time);
      tagms_adc_ped->Fill(pedestal_avg);
      tagms_adc_qf->Fill((*iter)->QF);
      tagms_adc_nsi->Fill((*iter)->nsamples_integral);
      tagms_adc_nsp->Fill((*iter)->nsamples_pedestal);

      // Fill 2D histograms
      tagms_adc_pint_2d->Fill(fiberNo, log10((*iter)->pulse_integral));
      tagms_adc_time_2d->Fill(fiberNo, (*iter)->pulse_time);
      tagms_adc_ped_2d->Fill(fiberNo, pedestal_avg);
      tagms_adc_qf_2d->Fill(fiberNo, (*iter)->QF);
      tagms_adc_nsi_2d->Fill(fiberNo, (*iter)->nsamples_integral);
      tagms_adc_nsp_2d->Fill(fiberNo, (*iter)->nsamples_pedestal);
      tagms_adc_pint_nsi->Fill((*iter)->nsamples_integral,
                              log10((*iter)->pulse_integral));
      tagms_adc_ped_nsp->Fill((*iter)->nsamples_pedestal,
                             pedestal_avg);

      // Fill dynamic array of histograms
      tagms_adc_pint_sng[fiberNo - 1]->Fill(log10((*iter)->pulse_integral));
      tagms_adc_time_sng[fiberNo - 1]->Fill((*iter)->pulse_time);
      tagms_adc_ped_sng[fiberNo - 1]->Fill(pedestal_avg);
      tagms_adc_qf_sng[fiberNo - 1]->Fill((*iter)->QF);

      // Calculate the FADC250 multiplicities for each channel
      ++single_adc_hits[fiberNo - 1];
      ++single_adc_hits_total;
    }
  }

  // DSCTDCDigiHits
  for (titer = tdcdigihits.begin(); titer != tdcdigihits.end(); ++titer) {
    int row = (*titer)->row;
    int column = (*titer)->column;
    int fiberNo = (row < 1)?      0 :
                  (column == 7)?  row :
                  (column == 25)? row + 5 :
                  (column == 79)? row + 10 :
                  (column == 97)? row + 15 : 0;
    if (fiberNo > 0) {

      // Fill 1D histograms
      tagms_tdc_seen->Fill(fiberNo);
      tagms_tdc_time->Fill((*titer)->time);

      // Fill 2D histograms
      tagms_tdc_time_2d->Fill(fiberNo, (*titer)->time);

      // Fill dynamic array of 1D histograms
      tagms_tdc_time_sng[fiberNo - 1]->Fill((*titer)->time);

      // Calculate the F1TDC multiplicities for each channel
      ++single_tdc_hits[fiberNo - 1];
      ++single_tdc_hits_total;
    }
  }

  // Fill the 2D multiplicity histograms for ADC and TDC
  for (unsigned int fno = 0; fno < NSINGLES; ++fno){
    if (single_adc_hits[fno] > 0) {
      tagms_adc_mult_2d->Fill(fno + 1, single_adc_hits[fno]);
      tagms_adc_mult_sng[fno]->Fill(single_adc_hits[fno]);
    }

    if (single_tdc_hits[fno] > 0) {
      tagms_tdc_mult_2d->Fill(fno + 1, single_tdc_hits[fno]);
      tagms_tdc_mult_sng[fno]->Fill(single_tdc_hits[fno]);
    }
  }

  // DTAGMHits
  for (hiter = hits.begin(); hiter != hits.end(); ++hiter) {
    int row = (*hiter)->row;
    int column = (*hiter)->column;
    int fiberNo = (row < 1)?      0 :
                  (column == 7)?  row :
                  (column == 25)? row + 5 :
                  (column == 79)? row + 10 :
                  (column == 97)? row + 15 : 0;
    if (fiberNo > 0) {

      // Fill 1D histograms
      tagms_hit_seen->Fill(fiberNo);
      tagms_hit_npix->Fill(log10((*hiter)->npix_fadc));
      tagms_hit_tadc->Fill((*hiter)->time_fadc);
      tagms_hit_time->Fill((*hiter)->t);
      tagms_hit_time_tadc->Fill((*hiter)->t, (*hiter)->time_fadc);
      tagms_hit_tdiff->Fill((*hiter)->t - (*hiter)->time_fadc);
 
      // Fill 2D histograms
      tagms_hit_npix_2d->Fill(fiberNo, (*hiter)->npix_fadc);
      tagms_hit_tadc_2d->Fill(fiberNo, (*hiter)->time_fadc);
      tagms_hit_time_2d->Fill(fiberNo, (*hiter)->t);

      // Fill dynamic array of 1D histograms
      tagms_hit_time_sng[fiberNo - 1]->Fill((*hiter)->t);
    }
  }

  // Find all possible time differences for multiple TDC hits
  if ((single_tdc_hits_total > 1) && (single_adc_hits_total > 1) ) {
    for (titer = tdcdigihits.begin(); titer != tdcdigihits.end(); ++titer) {
      for (iter = digihits.begin(); iter != digihits.end(); ++iter) {
        if ((*iter)->row > 0 && (*iter)->column == (*titer)->column &&
            (*iter)->row == (*titer)->row) 
        {
          tagms_tdc_tdiff_all->Fill((*titer)->time * TDC_NS_PER_COUNT -
                                   (*iter)->pulse_time * ADC_NS_PER_COUNT);
        }
      }
    }
  }

  // Lock ROOT mutex so other threads won't interfere 
  japp->RootUnLock();
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGM_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TAGM_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
