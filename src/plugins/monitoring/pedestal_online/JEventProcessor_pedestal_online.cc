// $Id$
//
//    File: JEventProcessor_pedestal_online.cc
// Created: Tue Feb 23 13:01:07 EST 2016
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#include <stdint.h>
#include <vector>
#include <cmath>

#include "JEventProcessor_pedestal_online.h"
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>

using namespace std;
using namespace jana;

#include <DAQ/DF1TDCHit.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/JEventSource_EVIO.h>
#include <TTAB/DTranslationTable.h>

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TROOT.h>


static const int highcratenum=100;
// root hist pointers
static TProfile *pedestal_vevent[highcratenum];
static TTree *pedestal_vtime_tree[highcratenum];
static TH1F *pedestal_vtime_hist[highcratenum];

int treetime, pednumsamps;
float pedmean, pedrms;

float treepedestal[highcratenum];

long int recentwalltime=0;

static const time_t periodlength=30;  // seconds
long int globalstarttime=0;
long int globalstoptime=0;
long int periodstarttime[highcratenum];



// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_pedestal_online());
  }
} // "C"


//------------------
// JEventProcessor_pedestal_online (Constructor)
//------------------
JEventProcessor_pedestal_online::JEventProcessor_pedestal_online()
{
	VERBOSE = 0;

	if(gPARMS){
		gPARMS->SetDefaultParameter("pedestal_online:VERBOSE", VERBOSE, "pedestal_online verbosity level");
	}
}

//------------------
// ~JEventProcessor_pedestal_online (Destructor)
//------------------
JEventProcessor_pedestal_online::~JEventProcessor_pedestal_online()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_pedestal_online::init(void)
{
	if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::init()\n");

	// create root folder for DAQ and cd to it, store main dir
	maindir = gDirectory;
	peddir = maindir->mkdir("pedestal");
	peddir->cd();

	// Initialise histograms and variables
	for (int i=0; i<highcratenum; i++) {
		periodstarttime[i] = 0;
		pedestal_vevent[i] = NULL;
		pedestal_vtime_tree[i] = NULL;
	}
	
	// back to main dir
	maindir->cd();

	return NOERROR;
}


//------------------
// brun
//------------------
jerror_t JEventProcessor_pedestal_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_pedestal_online::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// first look for an epics event
	vector<const DEPICSvalue*> depicsvalue;
	loop->Get(depicsvalue);
	if (depicsvalue.size()>0) {
		recentwalltime = (long int)depicsvalue[0]->timestamp;
		if (VERBOSE>=2) printf("JEventProcessor_pedestal_online::evnt %li found epics event at time %li\n",eventnumber,recentwalltime);
		// set the initial values
		if (globalstarttime==0) {
			globalstarttime=recentwalltime;
			for (int i=0; i<highcratenum; i++) periodstarttime[i] = globalstarttime;
			if (VERBOSE>=1) printf("\nsetting the global start time %li  %li\n\n",recentwalltime,globalstarttime);
		}
		return NOERROR;
	}

	vector<const Df250PulseData*> f250PDs;
	vector<const Df250PulseIntegral*> f250PIs;
	vector<const Df125PulseIntegral*> f125PIs;

	loop->Get(f250PDs);
	loop->Get(f250PIs);
	loop->Get(f125PIs);
	
	// Although we are only filling objects local to this plugin, the directory changes: Global ROOT lock
	japp->RootWriteLock();

	if (peddir!=NULL) peddir->cd();
	
	// Access F250 from Df250PulseIntegral object - pre-Fall 2016 data
	for(unsigned int i=0; i<f250PIs.size(); i++) {
		const Df250PulseIntegral *hit = f250PIs[i];
		int rocid = hit->rocid;

		if(rocid>=0 && rocid<100) {
			char cratename[255],title[255];
			// only use time if it is good
			if (recentwalltime>0) {
				if (pedestal_vtime_tree[rocid]==NULL) {
					if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating time tree and histogram for crate %i\n",rocid);
					sprintf(cratename,"pedestal_vtime_tree_crate%i",rocid);
					sprintf(title,"Crate %i avg pedestal (F250)",rocid);
					pedestal_vtime_tree[rocid] = new TTree(cratename,title);
					pedestal_vtime_tree[rocid]->Branch("time",&treetime);
					//pedestal_vtime_tree[rocid]->Branch("pedestal",&treepedestal[rocid]);
					pedestal_vtime_tree[rocid]->Branch("pedmean",&pedmean);
					pedestal_vtime_tree[rocid]->Branch("pedrms",&pedrms);
					pedestal_vtime_tree[rocid]->Branch("pednumsamps",&pednumsamps);
					sprintf(cratename,"pedestal_vtime_hist_crate%i",rocid);
					pedestal_vtime_hist[rocid] = new TH1F(cratename,title,100,0.,0.);
				}
				// keep a running stats of the pedestal
				if (hit->pedestal>0) pedestal_vtime_hist[rocid]->Fill(hit->pedestal);
				// if more than periodlength has elapsed, then end the average and fill the tree
				if (recentwalltime > periodstarttime[rocid]+periodlength) {
					if (VERBOSE>=3) printf("JEventProcessor_pedestal_online::evnt  filling tree crate %i\n",rocid);
					treetime = periodstarttime[rocid];
					pedmean = pedestal_vtime_hist[rocid]->GetMean();
					pedrms = pedestal_vtime_hist[rocid]->GetRMS();
					pednumsamps = pedestal_vtime_hist[rocid]->GetEntries();
                    if (VERBOSE>=3) printf("\t\t%li  %li  %li  %li  %8.4f  %8.4f  %7i\n",
										   periodstarttime[rocid], globalstoptime, recentwalltime,  periodlength, 
										   pedmean, pedrms, pednumsamps);
					pedestal_vtime_tree[rocid]->Fill();
					pedestal_vtime_hist[rocid]->Reset();
				}
				// advance to the next period
				while (recentwalltime > periodstarttime[rocid]+periodlength) {
					periodstarttime[rocid]+=periodlength;
					globalstoptime = periodstarttime[rocid];
				}
			}


			if (pedestal_vevent[rocid]==NULL) {
				if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating event histogram for crate %i\n",rocid);
				sprintf(cratename,"pedestal_vevent_crate%i",rocid);
				sprintf(title,"Crate %i avg pedestal (F250);event num;pedestal (all chan avg)",rocid);
				pedestal_vevent[rocid] = new TProfile(cratename,title,200,0.0,10000.0);
				pedestal_vevent[rocid]->SetStats(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
				pedestal_vevent[rocid]->SetCanExtend(TH1::kXaxis);
#else
				pedestal_vevent[rocid]->SetBit(TH1::kCanRebin);
#endif
			} 
			pedestal_vevent[rocid]->Fill(eventnumber,hit->pedestal);
		}
	}

	// Access F250 from Df250PulseData object - post-Fall 2016 data
	for(unsigned int i=0; i<f250PDs.size(); i++) {
		const Df250PulseData *hit = f250PDs[i];
		int rocid = hit->rocid;
		if(rocid>=0 && rocid<100) {
			char cratename[255],title[255];
			// only use time if it is good
			if (recentwalltime>0) {
				if (pedestal_vtime_tree[rocid]==NULL) {
					if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating time tree and histogram for crate %i\n",rocid);
					sprintf(cratename,"pedestal_vtime_tree_crate%i",rocid);
					sprintf(title,"Crate %i avg pedestal (F250)",rocid);
					pedestal_vtime_tree[rocid] = new TTree(cratename,title);
					pedestal_vtime_tree[rocid]->Branch("time",&treetime);
					//pedestal_vtime_tree[rocid]->Branch("pedestal",&treepedestal[rocid]);
					pedestal_vtime_tree[rocid]->Branch("pedmean",&pedmean);
					pedestal_vtime_tree[rocid]->Branch("pedrms",&pedrms);
					pedestal_vtime_tree[rocid]->Branch("pednumsamps",&pednumsamps);
					sprintf(cratename,"pedestal_vtime_hist_crate%i",rocid);
					pedestal_vtime_hist[rocid] = new TH1F(cratename,title,100,0.,0.);
				}
				// keep a running stats of the (single event) pedestal
				if (hit->pedestal>0) pedestal_vtime_hist[rocid]->Fill(hit->pedestal/hit->nsamples_pedestal);
				// if more than periodlength has elapsed, then end the average and fill the tree
				if (recentwalltime > periodstarttime[rocid]+periodlength) {
					if (VERBOSE>=3) printf("JEventProcessor_pedestal_online::evnt  filling tree crate %i\n",rocid);

					treetime = periodstarttime[rocid];
					pedmean = pedestal_vtime_hist[rocid]->GetMean();
					pedrms = pedestal_vtime_hist[rocid]->GetRMS();
					pednumsamps = pedestal_vtime_hist[rocid]->GetEntries();
                    if (VERBOSE>=3) printf("\t\t%li  %li  %li  %li  %8.4f  %8.4f  %7i\n",
										   periodstarttime[rocid], globalstoptime, recentwalltime,  periodlength, 
										   pedmean, pedrms, pednumsamps);
					pedestal_vtime_tree[rocid]->Fill();
					pedestal_vtime_hist[rocid]->Reset();
				}

				// advance to the next period
				while (recentwalltime > periodstarttime[rocid]+periodlength) {
					periodstarttime[rocid]+=periodlength;
					globalstoptime = periodstarttime[rocid];
				}
			}

			if (pedestal_vevent[rocid]==NULL) {
				if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating event histogram for crate %i\n",rocid);
				sprintf(cratename,"pedestal_vevent_crate%i",rocid);
				sprintf(title,"Crate %i avg pedestal (F250);event num;pedestal (all chan avg)",rocid);
				pedestal_vevent[rocid] = new TProfile(cratename,title,200,0.0,10000.0);
				pedestal_vevent[rocid]->SetStats(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
				pedestal_vevent[rocid]->SetCanExtend(TH1::kXaxis);
#else
				pedestal_vevent[rocid]->SetBit(TH1::kCanRebin);
#endif
			} 

			// monitor single-sample pedestal
			double locPedestalFraction = (hit->nsamples_pedestal == 0) ? 0.0 : double(hit->pedestal)/double(hit->nsamples_pedestal);
			pedestal_vevent[rocid]->Fill(eventnumber, locPedestalFraction);
			pedestal_vevent[rocid]->Fill(eventnumber,0);
		}
	}

	// Access F125 from Df125PulseIntegral object
	for(unsigned int i=0; i<f125PIs.size(); i++) {
		const Df125PulseIntegral *hit = f125PIs[i];
		int rocid = hit->rocid;

		if(rocid>=0 && rocid<100) {
			char cratename[255],title[255];
			// only use time if it is good
			if (recentwalltime>0) {
				if (pedestal_vtime_tree[rocid]==NULL) {
					if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating time tree and histogram for crate %i\n",rocid);
					sprintf(cratename,"pedestal_vtime_tree_crate%i",rocid);
					sprintf(title,"Crate %i avg pedestal (F250)",rocid);
					pedestal_vtime_tree[rocid] = new TTree(cratename,title);
					pedestal_vtime_tree[rocid]->Branch("time",&treetime);
					//pedestal_vtime_tree[rocid]->Branch("pedestal",&treepedestal[rocid]);
					pedestal_vtime_tree[rocid]->Branch("pedmean",&pedmean);
					pedestal_vtime_tree[rocid]->Branch("pedrms",&pedrms);
					pedestal_vtime_tree[rocid]->Branch("pednumsamps",&pednumsamps);
					sprintf(cratename,"pedestal_vtime_hist_crate%i",rocid);
					pedestal_vtime_hist[rocid] = new TH1F(cratename,title,100,0.,0.);
				}
				// keep a running stats of the pedestal
				if (hit->pedestal>0) pedestal_vtime_hist[rocid]->Fill(hit->pedestal);
				// if more than periodlength has elapsed, then end the average and fill the tree
				if (recentwalltime > periodstarttime[rocid]+periodlength) {
					if (VERBOSE>=3) printf("JEventProcessor_pedestal_online::evnt  filling tree crate %i\n",rocid);
					treetime = periodstarttime[rocid];
					pedmean = pedestal_vtime_hist[rocid]->GetMean();
					pedrms = pedestal_vtime_hist[rocid]->GetRMS();
					pednumsamps = pedestal_vtime_hist[rocid]->GetEntries();
					if (VERBOSE>=3) printf("\t\t%li  %li  %li  %li  %f  %f  %i\n",
										   periodstarttime[rocid], globalstoptime, recentwalltime,  periodlength, 
										   pedmean, pedrms, pednumsamps);
					pedestal_vtime_tree[rocid]->Fill();
					pedestal_vtime_hist[rocid]->Reset();
				}
				// advance to the next period
				while (recentwalltime > periodstarttime[rocid]+periodlength) {
					periodstarttime[rocid]+=periodlength;
					globalstoptime = periodstarttime[rocid];
				}
			}

			if (pedestal_vevent[rocid]==NULL) {
				if (VERBOSE>=1) printf("JEventProcessor_pedestal_online::evnt  creating event histogram for crate %i\n",rocid);
				sprintf(cratename,"pedestal_vevent_crate%i",rocid);
				sprintf(title,"Crate %i avg pedestal (F125);event num;pedestal (all chan avg)",rocid);
				pedestal_vevent[rocid] = new TProfile(cratename,title,200,0.0,10000.0);
				pedestal_vevent[rocid]->SetStats(0);
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
				pedestal_vevent[rocid]->SetCanExtend(TH1::kXaxis);
#else
				pedestal_vevent[rocid]->SetBit(TH1::kCanRebin);
#endif
			} 
			pedestal_vevent[rocid]->Fill(eventnumber,hit->pedestal);
		}
	}

	maindir->cd();
	// Unlock ROOT
	japp->RootUnLock();

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t JEventProcessor_pedestal_online::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.

	// Lock ROOT
	japp->RootWriteLock();

	peddir->cd();

	char cratename[255],title[255];

	int timeofrun = (globalstoptime - globalstarttime);
	int numberofbins = timeofrun/periodlength;

	for (int i=0; i<highcratenum; i++) {
		// if there's a tree make a histogram
		if (pedestal_vtime_tree[i]!=NULL) {
			sprintf(title,"Crate %i avg pedestal (F250)",i);
			sprintf(cratename,"pedestal_vtime_hist_crate%i",i);
			pedestal_vtime_hist[i] = new TH1F(cratename,title,numberofbins,globalstarttime,globalstoptime);
			pedestal_vtime_hist[i]->SetStats(0);
			pedestal_vtime_hist[i]->GetXaxis()->SetTimeDisplay(1);
			pedestal_vtime_hist[i]->GetXaxis()->SetTimeFormat("%H:%M %F 1970-01-01 00:00:00");
			//			pedestal_vtime_hist[i]->GetXaxis()->SetTimeFormat("%m/%d %H:%M %F 1970-01-01 00:00:00");

			Int_t nevent = pedestal_vtime_tree[i]->GetEntries();
			for (Int_t event=0;event<nevent;event++) {
				pedestal_vtime_tree[i]->GetEvent(event);
				if (pednumsamps>5) {
					int bin = pedestal_vtime_hist[i]->GetXaxis()->FindBin(treetime);
					float pederr = pedrms/sqrt(pednumsamps);
					if (VERBOSE>=3) printf("crate %i event %i bin %i   %f  %f  %i  %f\n",i, event,bin,pedmean, pedrms, pednumsamps,pederr);
					pedestal_vtime_hist[i]->SetBinContent(bin,pedmean);
					pedestal_vtime_hist[i]->SetBinError(bin,pederr);
				}
			}
			pedestal_vtime_hist[i]->SetMinimum(pedestal_vtime_hist[i]->GetMinimum(0.1));
		}

		if (pedestal_vevent[i] != NULL) {
			pedestal_vevent[i]->SetMinimum(pedestal_vevent[i]->GetMinimum(0.1));	
		}
	}

	// back to main dir
	maindir->cd();

	// Unlock ROOT
	japp->RootUnLock();
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_pedestal_online::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

