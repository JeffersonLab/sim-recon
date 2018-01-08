// Plugin for the tagger microscope time-walk corrections
// Author: aebarnes

#include "JEventProcessor_TAGM_TW.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ROOT header files
#include <TH1.h>
#include <TH2.h>
#include <TDirectory.h>
// GlueX header files
#include "TAGGER/DTAGMHit.h"
#include "RF/DRFTime_factory.h"

// Define constants
const uint32_t NCOLUMNS = 102;
const uint32_t NSINGLES = 20;
const uint32_t NROWS = 5;

const int32_t PMIN = 0;
const int32_t PMAX = 2000;
const int32_t PBIN = (PMAX - PMIN)/16;

const int32_t TMIN = -50;
const int32_t TMAX = 50;
const int32_t TBIN = (TMAX - TMIN)/0.05;

const double TMIN_RF = -2.0;
const double TMAX_RF = 2.0;
const double TBIN_RF = (TMAX_RF - TMIN_RF)/0.05;

const double TMIN_TW = -10.0;
const double TMAX_TW = 15.0;
const double TBIN_TW = (TMAX_TW - TMIN_TW)/0.05;

// Define histograms
//    timewalk
static TH2I* h_dt_vs_pp[NCOLUMNS];
static TH2I* h_dt_vs_pp_tdc[NCOLUMNS];
static TH2I* h_dt_vs_pp_ind[NROWS][4];
static TH2I* h_dt_vs_pp_tdc_ind[NROWS][4];
//    tdc - adc
static TH2I* h_tdc_adc_all;
static TH2I* h_tdc_adc_all_ind;
static TH2I* h_t_adc_all;
static TH2I* h_t_adc_all_ind;
//    adc - rf
static TH2I* h_adc_rf_all;
static TH2I* h_adc_rf_all_ind;

// Define RFTime_factory
DRFTime_factory* dRFTimeFactory;

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_TAGM_TW());
}
} // "C"


//------------------
// JEventProcessor_TAGM_TW (Constructor)
//------------------
JEventProcessor_TAGM_TW::JEventProcessor_TAGM_TW()
{

}

//------------------
// ~JEventProcessor_TAGM_TW (Destructor)
//------------------
JEventProcessor_TAGM_TW::~JEventProcessor_TAGM_TW()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TAGM_TW::init(void)
{
   TDirectory *main = gDirectory;

   TDirectory *tagmDir = gDirectory->mkdir("TAGM_TW");
   tagmDir->cd();

   gDirectory->mkdir("tdc-rf")->cd();
   for (uint32_t i = 0; i < NCOLUMNS; ++i)
   {
      h_dt_vs_pp_tdc[i] = new TH2I(Form("h_dt_vs_pp_tdc_%i", i+1),
                               Form("#delta t vs. pulse peak for TAGM column %i;\
                               Pulse peak (adc counts);TDC - RF (ns)", i+1), \
                               PBIN, PMIN, PMAX, TBIN, TMIN, TMAX);
   }
   for (uint32_t i = 0; i < NROWS; ++i)
   {
      for (uint32_t j = 0; j < 4; ++j)
      {
         h_dt_vs_pp_tdc_ind[i][j] = new TH2I(Form("h_dt_vs_pp_tdc_ind_%i_%i", i+1, j+1),
                                        Form("#delta t vs. pulse peak for TAGM ind. column %i, row %i;\
                                        Pulse peak (adc counts);TDC - RF (ns)", j+1, i+1),
                                        PBIN, PMIN, PMAX, TBIN, TMIN, TMAX);

      }
   }
   tagmDir->cd();

   gDirectory->mkdir("t-rf")->cd();
   for (uint32_t i = 0; i < NCOLUMNS; ++i)
   {
      h_dt_vs_pp[i] = new TH2I(Form("h_dt_vs_pp_%i", i+1),
                               Form("#delta t vs. pulse peak for TAGM column %i;\
                               Pulse peak (adc counts);T - RF (ns)", i+1),\
                               PBIN, PMIN, PMAX, TBIN_TW, TMIN_TW, TMAX_TW);
   }
   for (uint32_t i = 0; i < NROWS; ++i)
   {
      for (uint32_t j = 0; j < 4; ++j)
      {
         h_dt_vs_pp_ind[i][j] = new TH2I(Form("h_dt_vs_pp_ind_%i_%i", i+1, j+1),
                                    Form("#delta t vs. pulse peak for TAGM ind. column %i, row %i;\
                                    Pulse peak (adc counts);T - RF (ns)", j+1, i+1),
                                    PBIN, PMIN, PMAX, TBIN_TW, TMIN_TW, TMAX_TW);

      }
   }
   tagmDir->cd();

   h_tdc_adc_all =      new TH2I("tdc_adc_all","Summed channels;TDC - ADC (ns);Column",
                                  TBIN/10, TMIN, TMAX, NCOLUMNS,1, NCOLUMNS+1);
   h_tdc_adc_all_ind =  new TH2I("tdc_adc_all_ind","Individual channels;TDC - ADC (ns);Column",
                                  TBIN/10, TMIN, TMAX, NSINGLES,1, NSINGLES+1);
   h_t_adc_all =        new TH2I("t_adc_all","Summed channels;TDC - ADC (ns);Column",
                                  TBIN/10, TMIN, TMAX, NCOLUMNS,1, NCOLUMNS+1);
   h_t_adc_all_ind =    new TH2I("t_adc_all_ind","Individual channels;TDC - ADC (ns);Column",
                                  TBIN/10, TMIN, TMAX, NSINGLES,1, NSINGLES+1);
   h_adc_rf_all =       new TH2I("adc_rf_all","Summed channels;ADC - RF (ns);Column",
                                  TBIN_RF, TMIN_RF, TMAX_RF, NCOLUMNS,1, NCOLUMNS+1);
   h_adc_rf_all_ind =   new TH2I("adc_rf_all_ind","Individual channels;ADC - RF (ns);Column",
                                  TBIN_RF, TMIN_RF, TMAX_RF, NSINGLES,1, NSINGLES+1);

   main->cd();

   return NOERROR;
	
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TAGM_TW::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   // Initialize RF time factory
   dRFTimeFactory = static_cast<DRFTime_factory*>(eventLoop->GetFactory("DRFTime"));

   // be sure that DRFTime_factory::init() and brun() are called
   vector<const DRFTime*> locRFTimes;
   eventLoop->Get(locRFTimes);

   return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TAGM_TW::evnt(JEventLoop *loop, uint64_t eventnumber)
{
   // FILL HISTOGRAMS
   // Since we are filling histograms local to this plugin, it will not interfere
   // with other ROOT operations: can use plugin-wide ROOT fill lock

   vector<const DTAGMHit*>      hits;
   loop->Get(hits, "Calib");

   vector<const DRFTime*>	locRFTimes;
   loop->Get(locRFTimes, "TAGH");
   const DRFTime* locRFTime = NULL;

   if (locRFTimes.size() > 0)
      locRFTime = locRFTimes[0];
   else
      return NOERROR;

   japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

   for (uint32_t i = 0; i < hits.size(); ++i) {
      if (!hits[i]->has_TDC || !hits[i]->has_fADC) continue;
      int col      = hits[i]->column;
      int row      = hits[i]->row;
      double pp    = hits[i]->pulse_peak;
      double tm_t  = hits[i]->t;
      double adc_t = hits[i]->time_fadc;
      double tdc_t = hits[i]->time_tdc;
      double rf_adc = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime, adc_t);

      if (row == 0)
      {
         h_dt_vs_pp[col-1]->Fill(pp, tm_t - rf_adc);
         h_dt_vs_pp_tdc[col-1]->Fill(pp, tdc_t - rf_adc);
         h_tdc_adc_all->Fill(tdc_t - adc_t, col);
         h_t_adc_all->Fill(tm_t - adc_t, col);
         h_adc_rf_all->Fill(adc_t - rf_adc, col);
      }
      else
      {
         if (col == 9)
         {
            h_dt_vs_pp_ind[row-1][0]->Fill(pp,tm_t - rf_adc);
            h_dt_vs_pp_tdc_ind[row-1][0]->Fill(pp,tdc_t - rf_adc);
            h_tdc_adc_all_ind->Fill(tdc_t - adc_t, row);
            h_t_adc_all_ind->Fill(tm_t - adc_t, row);
            h_adc_rf_all_ind->Fill(adc_t - rf_adc, row);
         }
         else if (col == 27)
         {
            h_dt_vs_pp_ind[row-1][1]->Fill(pp,tm_t - rf_adc);
            h_dt_vs_pp_tdc_ind[row-1][1]->Fill(pp,tdc_t - rf_adc);
            h_tdc_adc_all_ind->Fill(tdc_t - adc_t, row + 5);
            h_t_adc_all_ind->Fill(tm_t - adc_t, row + 5);
            h_adc_rf_all_ind->Fill(adc_t - rf_adc, row + 5);
         }
         else if (col == 81)
         {
            h_dt_vs_pp_ind[row-1][2]->Fill(pp,tm_t - rf_adc);
            h_dt_vs_pp_tdc_ind[row-1][2]->Fill(pp,tdc_t - rf_adc);
            h_tdc_adc_all_ind->Fill(tdc_t - adc_t, row + 10);
            h_t_adc_all_ind->Fill(tm_t - adc_t, row + 10);
            h_adc_rf_all_ind->Fill(adc_t - rf_adc, row + 10);
         }
         else if (col == 99)
         {
            h_dt_vs_pp_ind[row-1][3]->Fill(pp,tm_t - rf_adc);
            h_dt_vs_pp_tdc_ind[row-1][3]->Fill(pp,tdc_t - rf_adc);
            h_tdc_adc_all_ind->Fill(tdc_t - adc_t, row + 15);
            h_t_adc_all_ind->Fill(tm_t - adc_t, row + 15);
            h_adc_rf_all_ind->Fill(adc_t - rf_adc, row + 15);
         }
      }
   }

   japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TAGM_TW::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TAGM_TW::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}
