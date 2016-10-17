// Plugin for the tagger microscope time-walk corrections
// Author: aebarnes

#include "JEventProcessor_TAGM_TW.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ROOT header files
#include <TH2.h>
// TAGM header files
#include "TAGGER/DTAGMHit.h"

// Define constants
const uint32_t NCOLUMNS = 102;
const uint32_t NROWS = 5;

// Define histograms
static TH2I* h_dt_vs_pp[NCOLUMNS];
static TH2I* h_dt_vs_pp_ind[NROWS][4];

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
   // Name histograms 
   for (uint32_t i = 0; i < NCOLUMNS; ++i)
   {
      h_dt_vs_pp[i] = new TH2I(Form("h_dt_vs_pp_%i",i+1),
                               Form("Time difference vs. pulse peak for TAGM column %i;\
                               Pulse peak (adc counts);TDC - ADC (ns)",i+1),\
                               1000,0,2000,250,-10,15);
   }
   for (uint32_t i = 0; i < NROWS; ++i)
   {
      for (uint32_t j = 0; j < 4; ++j)
      {
         h_dt_vs_pp_ind[i][j] = new TH2I(Form("h_dt_vs_pp_ind_%i_%i",i+1,j+1),
                                  Form("Time difference vs. pulse peak for TAGM ind. column %i, row %i;\
                                  Pulse peak (adc counts);TDC - ADC (ns)",j+1,i+1),
                                  1000,0,2000,250,-10,15);

      }
   }

   return NOERROR;
	
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TAGM_TW::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TAGM_TW::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   vector<const DTAGMHit*>      hits;
   loop->Get(hits, "Calib");

   // FILL HISTOGRAMS
   // Since we are filling histograms local to this plugin, it will not interfere
   // with other ROOT operations: can use plugin-wide ROOT fill lock
   japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

   for (uint32_t i = 0; i < hits.size(); ++i) {
      if (!hits[i]->has_TDC || !hits[i]->has_fADC) continue;
      int col      = hits[i]->column;
      int row      = hits[i]->row;
      double pp    = hits[i]->pulse_peak;
      double tm_t  = hits[i]->t;
      double adc_t = hits[i]->time_fadc;

      double tdiff = tm_t - adc_t;

      if (row == 0)
         h_dt_vs_pp[col-1]->Fill(pp,tdiff);
      else
      {
         if (col == 9)
            h_dt_vs_pp_ind[row-1][0]->Fill(pp,tdiff);
         else if (col == 27)
            h_dt_vs_pp_ind[row-1][1]->Fill(pp,tdiff);
         else if (col == 81)
            h_dt_vs_pp_ind[row-1][2]->Fill(pp,tdiff);
         else if (col == 99)
            h_dt_vs_pp_ind[row-1][3]->Fill(pp,tdiff);
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
