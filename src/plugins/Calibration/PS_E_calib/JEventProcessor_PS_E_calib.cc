// $Id$
// //    File: JEventProcessor_PS_E_calib.cc // Created: Thu Jul  9 17:44:32 EDT 5015
// Creator: aebarnes (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_PS_E_calib.h"
using namespace jana;

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TProfile.h>
#include <TH2.h>

#include <TAGGER/DTAGMHit.h>
#include <TAGGER/DTAGHHit.h>
#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <PAIR_SPECTROMETER/DPSPair.h>

#define CORRECTIONS false

// Define constants
const float Ebw_PS = 0.013;			// Energy bin width for the PS in GeV
const float Ebl_PS = 2.3;			// Low energy of the PS total energy in GeV
const float Ebh_PS = 4.9;			// High energy of the PS total energy in GeV
const float NEb_PS = (Ebh_PS - Ebl_PS)/Ebw_PS;	// Number of energy bins for the PS

const int MAX_COLUMNS = 100;			// Total columns in the TAGM
const int MAX_COUNTERS = 274;			// Total possible counters in the TAGH

// Declare variables
double p0 = 0;					// PS energy correction parameter
double p1 = 0;					// PS energy correction parameter
double p2 = 0;					// PS energy correction parameter
// TAGM
int column = 0;					// TAGM column
double tm_E = 0;				// TAGM energy
double tm_t = 0;				// TAGM time
double tdiff_tm = 0;				// time difference PS avg - TAGM
// TAGH
int counter = 0;				// TAGH counter
double th_E = 0;				// TAGH energy
double th_t = 0;				// TAGH time
double tdiff_th = 0;				// time difference PS avg - TAGH
// PS
double ps_E = 0;				// total PS energy
double ps_El = 0;				// left PS arm energy
double ps_Er = 0;				// right PS arm energy
double ps_t = 0;				// average PS time
double ps_tl = 0;				// left PS time
double ps_tr = 0;				// right PS time

// This is not correct. Fix this when the timing is aligned
int run = 0;					// Run number

// Declare histograms
// TAGM
static TH2F *h_psE_vs_psEl_tm[MAX_COLUMNS];	// PS total vs fraction of PS left
// TAGH
static TH2F *h_psE_vs_psEl_th[MAX_COUNTERS];	// PS total vs fraction of PS left
// Timing check
static TH1F *h_dt;				// PS - TAGX timing check

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_PS_E_calib());
}
} // "C"


//------------------
// JEventProcessor_PS_E_calib (Constructor)
//------------------
JEventProcessor_PS_E_calib::JEventProcessor_PS_E_calib()
{

}

//------------------
// ~JEventProcessor_PS_E_calib (Destructor)
//------------------
JEventProcessor_PS_E_calib::~JEventProcessor_PS_E_calib()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_PS_E_calib::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

   // create root folder tagm
   TDirectory *mainDir = gDirectory;
   TDirectory *tagmDir = gDirectory->mkdir("TAGM");
   TDirectory *taghDir = gDirectory->mkdir("TAGH");

   // Book 2d histograms
   // TAGM
   tagmDir->cd();

   h_dt = new TH1F("h_dt","Time difference PS - TAGM;PS-TAGM (ns)",200,-17,37);

   for (int col = 0; col < MAX_COLUMNS; ++col) {
      h_psE_vs_psEl_tm[col] = new TH2F(Form("h_psE_vs_psEl_tm_%i",col+1),
                                       Form("PS E vs PS left, TAGM col %i;\
                                       Energy asymmetry;PS energy (GeV)",col+1),
                                       50,0,1,NEb_PS,Ebl_PS,Ebh_PS);
   }

   // TAGH
   taghDir->cd();

   for (int hodo = 0; hodo < MAX_COUNTERS; ++hodo) {
      h_psE_vs_psEl_th[hodo] = new TH2F(Form("psE_vs_psEl_th_%i",hodo+1),
                                        Form("PS E vs PS left, TAGH counter %i;\
                                        Energy asymmetry;PS energy (GeV)",hodo+1),
                                        50,0,1,NEb_PS,Ebl_PS,Ebh_PS);
   }
 

   mainDir->cd();
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_PS_E_calib::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes

   // THIS IS TEMPORARY! FIX ONCE PS-TAGX COINCIDENCE
   // IS THE SAME FOR ALL RUNS
   run = runnumber;

   // Get the PS energy corrections from CCDB 
   std::vector< std::map<std::string, double> > table;
   std::string ccdb_key = "/PHOTON_BEAM/pair_spectrometer/fine/energy_corrections";
   if (eventLoop->GetCalib(ccdb_key, table)) {
      jout << "Error loading " << ccdb_key << " from ccdb!" << std::endl;
   }
   for (unsigned int i=0; i < table.size(); ++i) {
      p0 = (table[i])["constant"];
      p1 = (table[i])["linear"];
      p2 = (table[i])["quadratic"];
   }

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_PS_E_calib::evnt(JEventLoop *loop, uint64_t eventnumber)
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

   vector<const DTAGMHit*>	tm_hits;
   vector<const DTAGHHit*>	th_hits;
   vector<const DPSPair*>	ps_pairs;
   vector<const DPSCPair*>	psc_pairs;

   loop->Get(tm_hits);
   loop->Get(th_hits);
   loop->Get(ps_pairs);
   loop->Get(psc_pairs);

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

   if (psc_pairs.size() > 0 ) {
      for (uint32_t i = 0; i < ps_pairs.size(); ++i) {
         // left arm
         ps_El = ps_pairs[i]->ee.first->E;
         ps_tl = ps_pairs[i]->ee.first->t;
         // right arm
         ps_Er = ps_pairs[i]->ee.second->E;
         ps_tr = ps_pairs[i]->ee.second->t;
         // combined
         ps_t = (ps_tl + ps_tr)/2;
         ps_E = ps_El + ps_Er;

         #if CORRECTIONS
         ps_E /= p0 + p1*(ps_El/ps_E) + p2*(ps_El/ps_E)*(ps_El/ps_E);
         #endif
         
         // loop over TAGM hits
         for (uint32_t j = 0; j < tm_hits.size(); ++j) {
            tm_E = tm_hits[j]->E;
            tm_t = tm_hits[j]->t;
            tdiff_tm = ps_t - tm_t;
            column = tm_hits[j]->column;

            if (!tm_hits[j]->has_fADC || !tm_hits[j]->has_TDC) continue;
            if (tm_hits[j]->row != 0) continue;

            h_dt->Fill(tdiff_tm);

            if (run == 3180 && tdiff_tm < 26.5 && tdiff_tm >= 22.5) {
               h_psE_vs_psEl_tm[column-1]->Fill(ps_El/ps_E,ps_E);
            }
            else if (run == 3185 && tdiff_tm < 34.5 && tdiff_tm >= 30.5) {
               h_psE_vs_psEl_tm[column-1]->Fill(ps_El/ps_E,ps_E);
            }
         }
         
         // loop over TAGH hits
         for (uint32_t j = 0; j < th_hits.size(); ++j) {
            th_E = th_hits[j]->E;
            th_t = th_hits[j]->t;
            tdiff_th = ps_t - th_t;
            counter = th_hits[j]->counter_id;

            if (!th_hits[j]->has_fADC || !th_hits[j]->has_TDC) continue;

            if (run == 3180 && tdiff_th < 26.5 && tdiff_th >= 22.5) {
               h_psE_vs_psEl_th[counter-1]->Fill(ps_El/ps_E,ps_E);
            }
            else if (run == 3185 && tdiff_th < 34.5 && tdiff_th >= 30.5) {
               h_psE_vs_psEl_th[counter-1]->Fill(ps_El/ps_E,ps_E);
            }
         }
      }
   }

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_PS_E_calib::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_PS_E_calib::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

