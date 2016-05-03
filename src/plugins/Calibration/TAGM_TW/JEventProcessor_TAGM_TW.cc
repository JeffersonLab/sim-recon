// Plugin for the tagger microscope time-walk corrections
// Author: aebarnes

#include "JEventProcessor_TAGM_TW.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ROOT header files
#include <TH2.h>
// C++ header files
#include <stdint.h>
#include <vector>
#include <stdio.h>
// TAGM header files
#include "TAGGER/DTAGMHit.h"
// PS header files
#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <PAIR_SPECTROMETER/DPSPair.h>
// RF header files
#include <RF/DRFTime_factory.h>

// Define constants
const uint32_t NCOLUMNS = 100;
Double_t p0 = 0;	// Parameter for PS energy calibration
Double_t p1 = 0;	// Parameter for PS energy calibration
Double_t p2 = 0;	// Parameter for PS energy calibration
const Double_t c0 = 3.52688;	// Parameter for PS energy conversion to TAGH scale
const Double_t c1 = -1.11307;	// Parameter for PS energy conversion to TAGH scale
const Double_t c2 = 0.319561;	// Parameter for PS energy conversion to TAGH scale

// Define histograms
static TH2F* h_dt_vs_pp[NCOLUMNS];

// Define variables
Double_t tm_t;		// TAGM tdc time
Double_t psc_t;		// psc tdc time
Double_t rf_t;		// rf time from the TAGH RF source
Double_t tdiff_rf;	// Time difference TAGM - RF
Double_t tm_E;		// TAGM energy
Double_t ps_E;		// PS energy
Double_t ps_El;		// PS left arm energy
Double_t ps_Er;		// PS right arm energy
Double_t adc_pp; 	// pulse peak
Int_t col;		// TAGM column
Int_t row;		// TAGM row

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
   // Name histograms 
   for (uint32_t i = 0; i < NCOLUMNS; ++i) {
      h_dt_vs_pp[i] = new TH2F(Form("h_dt_vs_pp_%i",i+1),
                               Form("Time difference vs. pulse peak for TAGM column %i;\
                               Pulse peak (adc counts);TAGM - RF (ns)",i+1),\
                               1000,0,1000,200,-37,-17);
   }

   return NOERROR;
	
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TAGM_TW::brun(JEventLoop *eventLoop, int32_t runnumber)
{
   //////////
   // RF
   // ///////

   // Initialize RF time factory
   dRFTimeFactory = static_cast<DRFTime_factory*>(eventLoop->GetFactory("DRFTime"));

   // be sure that DRFTime_factory::init() and brun() are called
   vector<const DRFTime*> locRFTimes;

   eventLoop->Get(locRFTimes);

   // Load PS energy correction values from CCDB
   std::vector< std::map<std::string, double> > table;
   std::string ccdb_key = "/PHOTON_BEAM/pair_spectrometer/fine/energy_corrections";
   if (eventLoop->GetCalib(ccdb_key, table)) {
      jout << "Error loading " << ccdb_key << " from ccdb!" << std::endl;
   }
   p0 = (table[0])["constant"];
   p1 = (table[0])["linear"];
   p2 = (table[0])["quadratic"];

   return NOERROR;

}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TAGM_TW::evnt(JEventLoop *loop, uint64_t eventnumber)
{

   vector<const DTAGMHit*>      hits;
   vector<const DPSCPair*>	psc_pairs;
   vector<const DPSPair*>	ps_pairs;
   vector<const DRFTime*>	locRFTimes;

   loop->Get(hits);
   loop->Get(psc_pairs);
   loop->Get(ps_pairs);
   loop->Get(locRFTimes, "PSC");
   const DRFTime* locRFTime = NULL;

   if (locRFTimes.size() > 0)
      locRFTime = locRFTimes[0];
   else
      return NOERROR;

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

   if (psc_pairs.size() > 0 && ps_pairs.size() > 0) {
      for (uint32_t i = 0; i < hits.size(); ++i) {
         if (!hits[i]->has_TDC || !hits[i]->has_fADC) continue;
         // TAGM
         col    = hits[i]->column;
         row    = hits[i]->row;
         adc_pp = hits[i]->pulse_peak;
         tm_t   = hits[i]->t;
         tm_E   = hits[i]->E;

         // PS
         ps_El = ps_pairs[0]->ee.first->E;
         ps_Er = ps_pairs[0]->ee.second->E;
         ps_E  = ps_El + ps_Er;

         // Apply PS energy calibration
         ps_E /= p0 + p1*(ps_El/ps_E) + p2*(ps_El/ps_E)*(ps_El/ps_E);

         // Apply PS conversion to TAGH energy scale
         ps_E = c0 + c1*ps_E + c2*ps_E*ps_E;

         // Use PS time to get RF time
         psc_t  = psc_pairs[0]->ee.first->t;
         rf_t   = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime, psc_t);
         tdiff_rf = tm_t - rf_t;

         if (tm_E <= ps_E*1.01 && tm_E >= ps_E*0.99 && row == 0) {
            h_dt_vs_pp[col - 1]->Fill(adc_pp,tdiff_rf);
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
