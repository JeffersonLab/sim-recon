// $Id$
//
//    File: JEventProcessor_TAGM_thresh.cc
// Created: Tue Feb 28 17:43:55 EST 2017
// Creator: barnes (on Linux gluey.phys.uconn.edu 2.6.32-642.13.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_TAGM_thresh.h"
using namespace jana;

#include <TAGGER/DTAGMHit.h>
#include <TH1.h>
#include <TDirectory.h>

static const uint32_t NCOLUMNS = 102;
static const uint32_t NSINGLES = 20;
static const uint32_t NROWS = 5;

static const uint32_t IMIN = 0;
static const uint32_t IMAX = 5;
static const uint32_t IBINS = (IMAX-IMIN)/0.01;

static const uint32_t PMIN = 0;
static const uint32_t PMAX = 4095;
static const uint32_t PBINS = (PMAX-PMIN)/10;

static TH1I *h_int[NCOLUMNS];
static TH1I *h_int_ind[NSINGLES];
static TH1I *h_pp[NCOLUMNS];
static TH1I *h_pp_ind[NSINGLES];

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_TAGM_thresh());
}
} // "C"


//------------------
// JEventProcessor_TAGM_thresh (Constructor)
//------------------
JEventProcessor_TAGM_thresh::JEventProcessor_TAGM_thresh()
{

}

//------------------
// ~JEventProcessor_TAGM_thresh (Destructor)
//------------------
JEventProcessor_TAGM_thresh::~JEventProcessor_TAGM_thresh()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TAGM_thresh::init(void)
{
	// This is called once at program startup. 

   TDirectory *tagmDir = gDirectory->mkdir("TAGM_thresh");
   gDirectory->mkdir("TAGM_thresh/integrals");
   gDirectory->mkdir("TAGM_thresh/peaks");

   tagmDir->cd();
   for (uint32_t i = 0; i < NCOLUMNS; ++i)
   {
      tagmDir->cd("integrals");
      h_int[i] = new TH1I(Form("h_int_%i",i+1),
                          Form("Pulse integral for col %i;\
                                log10(pulse integral)",i+1),
                          IBINS,IMIN,IMAX);
      tagmDir->cd("peaks");
      h_pp[i] = new TH1I(Form("h_pp_%i",i+1),
                         Form("Pulse height for col %i;\
                               ADC",i+1),PBINS,PMIN,PMAX);
   }

   for (uint32_t i = 0; i < NSINGLES; ++i)
   {
      tagmDir->cd("integrals");
      h_int_ind[i] = new TH1I(Form("h_int_ind_%i",i+1),
                              Form("Pulse integral for ind. col %i;\
                                    log10(pulse integral)",i+1),
                              IBINS,IMIN,IMAX);
      tagmDir->cd("peaks");
      h_pp_ind[i] = new TH1I(Form("h_pp_ind_%i",i+1),
                             Form("Pulse height for ind. col %i;\
                                   ADC",i+1),PBINS,PMIN,PMAX);
   }
   tagmDir->cd();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TAGM_thresh::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TAGM_thresh::evnt(JEventLoop *loop, uint64_t eventnumber)
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
	// japp->RootFillLock(this);
	//  ... fill historgrams or trees ...
	// japp->RootFillUnLock(this);

   vector<const DTAGMHit*>	hits;
   loop->Get(hits, "Calib");

   japp->RootFillLock(this);
   for (uint32_t i = 0; i < hits.size(); ++i)
   {
      if (!hits[i]->has_fADC) continue;
      int col = hits[i]->column;
      int row = hits[i]->row;
      double pp = hits[i]->pulse_peak;
      double pi = hits[i]->integral;

      if (row == 0)
      {
         h_int[col-1]->Fill(log10(pi));
         h_pp[col-1]->Fill(pp);
      }
      else if (col == 9)
      {
         h_int_ind[row-1]->Fill(log10(pi));
         h_pp_ind[row-1]->Fill(pp);
      }
      else if (col == 27)
      {
         h_int_ind[row-1+5]->Fill(log10(pi));
         h_pp_ind[row-1+5]->Fill(pp);
      }
      else if (col == 81)
      {
         h_int_ind[row-1+10]->Fill(log10(pi));
         h_pp_ind[row-1+10]->Fill(pp);
      }
      else if (col == 99)
      {
         h_int_ind[row-1+15]->Fill(log10(pi));
         h_pp_ind[row-1+15]->Fill(pp);
      }
   }
   japp->RootFillUnLock(this);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TAGM_thresh::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TAGM_thresh::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

