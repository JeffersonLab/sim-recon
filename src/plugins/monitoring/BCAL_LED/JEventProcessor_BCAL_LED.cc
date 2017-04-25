// $Id$
//
//    File: JEventProcessor_BCAL_LED.cc
//

#include <stdint.h>
#include <vector>
#include "TTree.h"
#include "JEventProcessor_BCAL_LED.h"
#include <JANA/JApplication.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace jana;

#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250WindowRawData.h"
#include "TRIGGER/DL1Trigger.h"

#include <TDirectory.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TStyle.h>

// root hist pointers

//all channels
static TProfile *bcal_peak_vevent = NULL;

//2 sides
static TProfile *up_peak_vevent = NULL;
static TProfile *down_peak_vevent = NULL;
//4 columns
static TProfile *column1_peak_vevent = NULL;
static TProfile *column2_peak_vevent = NULL;
static TProfile *column3_peak_vevent = NULL;
static TProfile *column4_peak_vevent = NULL;

static TProfile *column1_up_peak_vevent = NULL;
static TProfile *column2_up_peak_vevent = NULL;
static TProfile *column3_up_peak_vevent = NULL;
static TProfile *column4_up_peak_vevent = NULL;
static TProfile *column1_down_peak_vevent = NULL;
static TProfile *column2_down_peak_vevent = NULL;
static TProfile *column3_down_peak_vevent = NULL;
static TProfile *column4_down_peak_vevent = NULL;


static TProfile *column1_up_peak_vevent1 = NULL;
static TProfile *column1_down_peak_vevent1 = NULL;
static TProfile *column1_up_peak_vevent2 = NULL;
static TProfile *column1_down_peak_vevent2 = NULL;
static TProfile *column1_up_peak_vevent3 = NULL;
static TProfile *column1_down_peak_vevent3 = NULL;
static TProfile *column1_up_peak_vevent4 = NULL;
static TProfile *column1_down_peak_vevent4 = NULL;

static TProfile *column2_up_peak_vevent1 = NULL;
static TProfile *column2_down_peak_vevent1 = NULL;
static TProfile *column2_up_peak_vevent2 = NULL;
static TProfile *column2_down_peak_vevent2 = NULL;
static TProfile *column2_up_peak_vevent3 = NULL;
static TProfile *column2_down_peak_vevent3 = NULL;
static TProfile *column2_up_peak_vevent4 = NULL;
static TProfile *column2_down_peak_vevent4 = NULL;

static TProfile *column3_up_peak_vevent1 = NULL;
static TProfile *column3_down_peak_vevent1 = NULL;
static TProfile *column3_up_peak_vevent2 = NULL;
static TProfile *column3_down_peak_vevent2 = NULL;
static TProfile *column3_up_peak_vevent3 = NULL;
static TProfile *column3_down_peak_vevent3 = NULL;
static TProfile *column3_up_peak_vevent4 = NULL;
static TProfile *column3_down_peak_vevent4 = NULL;

static TProfile *column4_up_peak_vevent1 = NULL;
static TProfile *column4_down_peak_vevent1 = NULL;
static TProfile *column4_up_peak_vevent2 = NULL;
static TProfile *column4_down_peak_vevent2 = NULL;
static TProfile *column4_up_peak_vevent3 = NULL;
static TProfile *column4_down_peak_vevent3 = NULL;
static TProfile *column4_up_peak_vevent4 = NULL;
static TProfile *column4_down_peak_vevent4 = NULL;

  
    
static TProfile *low_up_1 = NULL;
static TProfile *low_up_2 = NULL;
static TProfile *low_up_3 = NULL;
static TProfile *low_up_4 = NULL;
static TProfile *low_down_1 = NULL;
static TProfile *low_down_2 = NULL;
static TProfile *low_down_3 = NULL;
static TProfile *low_down_4 = NULL;

static TProfile *high_up_1 = NULL;
static TProfile *high_up_2 = NULL;
static TProfile *high_up_3 = NULL;
static TProfile *high_up_4 = NULL;
static TProfile *high_down_1 = NULL;
static TProfile *high_down_2 = NULL;
static TProfile *high_down_3 = NULL;
static TProfile *high_down_4 = NULL;

// Histograms added by Elton for z distributions



static TProfile* h2_ledboth_Aall_vs_event = NULL;
static TProfile* h2_ledboth_sector_vs_event = NULL;

static TH1I* h1_ledup_z_all = NULL;;
static TH2I* h2_ledup_z_vs_cellid = NULL;
static TH1I* h1_ledup_sector = NULL;
static TH1I* h1_ledup_sector_config = NULL;
static TH1I* h1_ledup_Tdiff_all = NULL;
static TH1I* h1_ledup_Tup_all = NULL;
static TH1I* h1_ledup_Tdown_all = NULL;
static TH1I* h1_ledup_Aup_all = NULL;
static TH1I* h1_ledup_Adown_all = NULL;
static TH2I* h2_ledup_Aup_vs_z = NULL;
static TH2I* h2_ledup_Adown_vs_z = NULL;
static TProfile* h2_ledup_Aup_vs_event = NULL;
static TProfile* h2_ledup_Adown_vs_event = NULL;
static TProfile* h2_ledup_Aall_vs_event = NULL;
static TProfile* h2_ledup_sector_vs_event = NULL;

static TH1I* h1_leddown_z_all = NULL;
static TH2I* h2_leddown_z_vs_cellid = NULL;
static TH1I* h1_leddown_sector = NULL;
static TH1I* h1_leddown_sector_config = NULL;
static TH1I* h1_leddown_Tdiff_all = NULL;
static TH1I* h1_leddown_Tup_all = NULL;
static TH1I* h1_leddown_Tdown_all = NULL;
static TH1I* h1_leddown_Aup_all = NULL;
static TH1I* h1_leddown_Adown_all = NULL;
static TH2I* h2_leddown_Aup_vs_z = NULL;
static TH2I* h2_leddown_Adown_vs_z = NULL;
static TProfile* h2_leddown_Aup_vs_event = NULL;
static TProfile* h2_leddown_Adown_vs_event = NULL;
static TProfile* h2_leddown_Aall_vs_event = NULL;
static TProfile* h2_leddown_sector_vs_event = NULL;




//----------------------------------------------------------------------------------

// 	string make_filename( const string& basename, int index, const string& ext )
// 	  {
// 	  ostringstream result;
// 	  result << basename << index << ext;
// 	  return result.str();
// 	  }


// Routine used to create our JEventProcessor
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_BCAL_LED());
	}

}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LED::JEventProcessor_BCAL_LED() {
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LED::~JEventProcessor_BCAL_LED() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_BCAL_LED::init(void) {
	
	// lock all root operations
	japp->RootWriteLock();
	
	// First thread to get here makes all histograms. If one pointer is
	// already not NULL, assume all histograms are defined and return now
	if(bcal_peak_vevent != NULL){
		japp->RootUnLock();
		return NOERROR;
	}
	
	//NOtrig=0; FPtrig=0; GTPtrig=0; FPGTPtrig=0; trigUS=0; trigDS=0; trigCosmic=0;
	//low_down_1_counter=0; low_down_2_counter=0; low_down_3_counter=0; low_down_4_counter=0; low_up_1_counter=0; low_up_2_counter=0; low_up_3_counter=0; 		low_up_4_counter=0; high_down_1_counter=0; high_down_2_counter=0; high_down_3_counter=0; high_down_4_counter=0; high_up_1_counter=0;
	//high_up_2_counter=0; high_up_3_counter=0; high_up_4_counter=0;
	//unidentified = 0; ledcounter = 0;

	adccount1 = 1100;
	adccount2 = 1200;
	
	
	maxnumberofevents=700000000.0;//Assuming 1Hz LED trigger, 300M for a beam run with 30KHz trigger and 700M for 70KHz
	//maxnumberofevents=10000.0;//using LED event conter
	nbins=24002;//Assuming 1Hz LED trigger, 10K for a beam run with 30KHz trigger and 24K for 70KHz
	//nbins=375002;//Based on cosmic run with 800Hz trigger and 1Hz LED trigger
	//nbins=3750002;//Based on cosmic run with 800Hz trigger and 10Hz LED trigger

	
	//overflow=0; underflow=0; negatives=0; zeros=0;

	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcalLED")->cd();
	//gStyle->SetOptStat(111110);


	
	bcal_peak_vevent = new TProfile("bcal_peak_vevent","Avg BCAL peak vs event;event num;peak (all chan avg)",nbins,0.0,maxnumberofevents);
	
	up_peak_vevent = new TProfile("up_peak_vevent","Avg BCAL peak vs event;event num;peak (all upstream chan avg)",nbins,0.0,maxnumberofevents);
	down_peak_vevent = new TProfile("down_peak_vevent","Avg BCAL peak vs event;event num;peak (all downstream chan avg)",nbins,0.0,maxnumberofevents);
	
	column1_peak_vevent = new TProfile("column1_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 1 chan avg)",nbins,0.0,maxnumberofevents);
	column2_peak_vevent = new TProfile("column2_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 2 chan avg)",nbins,0.0,maxnumberofevents);
	column3_peak_vevent = new TProfile("column3_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 3 chan avg)",nbins,0.0,maxnumberofevents);
	column4_peak_vevent = new TProfile("column4_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 4 chan avg)",nbins,0.0,maxnumberofevents);
	
	column1_up_peak_vevent = new TProfile("column1_up_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_up_peak_vevent = new TProfile("column2_up_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_up_peak_vevent = new TProfile("column3_up_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_up_peak_vevent = new TProfile("column4_up_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_down_peak_vevent = new TProfile("column1_down_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 1 down chan avg)",nbins,0.0,maxnumberofevents);
	column2_down_peak_vevent = new TProfile("column2_down_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 2 down chan avg)",nbins,0.0,maxnumberofevents);
	column3_down_peak_vevent = new TProfile("column3_down_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 3 down chan avg)",nbins,0.0,maxnumberofevents);
	column4_down_peak_vevent = new TProfile("column4_down_peak_vevent","Avg BCAL peak vs event;event num;peak (all column 4 down chan avg)",nbins,0.0,maxnumberofevents);
	
	column1_up_peak_vevent1 = new TProfile("column1_up_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_down_peak_vevent1 = new TProfile("column1_down_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 1 down chan avg)",nbins,0.0,maxnumberofevents);
	column1_up_peak_vevent2 = new TProfile("column1_up_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_down_peak_vevent2 = new TProfile("column1_down_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_up_peak_vevent3 = new TProfile("column1_up_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_down_peak_vevent3 = new TProfile("column1_down_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_up_peak_vevent4 = new TProfile("column1_up_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);
	column1_down_peak_vevent4 = new TProfile("column1_down_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 1 up chan avg)",nbins,0.0,maxnumberofevents);

	column2_up_peak_vevent1 = new TProfile("column2_up_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_down_peak_vevent1 = new TProfile("column2_down_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 2 down chan avg)",nbins,0.0,maxnumberofevents);
	column2_up_peak_vevent2 = new TProfile("column2_up_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_down_peak_vevent2 = new TProfile("column2_down_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_up_peak_vevent3 = new TProfile("column2_up_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_down_peak_vevent3 = new TProfile("column2_down_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_up_peak_vevent4 = new TProfile("column2_up_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);
	column2_down_peak_vevent4 = new TProfile("column2_down_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 2 up chan avg)",nbins,0.0,maxnumberofevents);

	column3_up_peak_vevent1 = new TProfile("column3_up_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_down_peak_vevent1 = new TProfile("column3_down_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 3 down chan avg)",nbins,0.0,maxnumberofevents);
	column3_up_peak_vevent2 = new TProfile("column3_up_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_down_peak_vevent2 = new TProfile("column3_down_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_up_peak_vevent3 = new TProfile("column3_up_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_down_peak_vevent3 = new TProfile("column3_down_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_up_peak_vevent4 = new TProfile("column3_up_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);
	column3_down_peak_vevent4 = new TProfile("column3_down_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 3 up chan avg)",nbins,0.0,maxnumberofevents);

	column4_up_peak_vevent1 = new TProfile("column4_up_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_down_peak_vevent1 = new TProfile("column4_down_peak_vevent1","Avg BCAL peak vs event;event num;peak (all column 4 down chan avg)",nbins,0.0,maxnumberofevents);
	column4_up_peak_vevent2 = new TProfile("column4_up_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_down_peak_vevent2 = new TProfile("column4_down_peak_vevent2","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_up_peak_vevent3 = new TProfile("column4_up_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_down_peak_vevent3 = new TProfile("column4_down_peak_vevent3","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_up_peak_vevent4 = new TProfile("column4_up_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);
	column4_down_peak_vevent4 = new TProfile("column4_down_peak_vevent4","Avg BCAL peak vs event;event num;peak (all column 4 up chan avg)",nbins,0.0,maxnumberofevents);

	
	
	low_up_1 = new TProfile("low_bias_up_column_1_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-4,1536);
	low_up_2 = new TProfile("low_bias_up_column_2_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-3,1537);
	low_up_3 = new TProfile("low_bias_up_column_3_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-2,1538);
	low_up_4 = new TProfile("low_bias_up_column_4_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-1,1539);
	
	low_down_1 = new TProfile("low_bias_down_column_1_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-4,1536);
	low_down_2 = new TProfile("low_bias_down_column_2_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-3,1537);
	low_down_3 = new TProfile("low_bias_down_column_3_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-2,1538);
	low_down_4 = new TProfile("low_bias_down_column_4_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-1,1539);
	
	high_up_1 = new TProfile("high_bias_up_column_1_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-4,1536);
	high_up_2 = new TProfile("high_bias_up_column_2_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-3,1537);
	high_up_3 = new TProfile("high_bias_up_column_3_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-2,1538);
	high_up_4 = new TProfile("high_bias_up_column_4_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-1,1539);
	
	high_down_1 = new TProfile("high_bias_down_column_1_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-4,1536);
	high_down_2 = new TProfile("high_bias_down_column_2_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-3,1537);
	high_down_3 = new TProfile("high_bias_down_column_3_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-2,1538);
	high_down_4 = new TProfile("high_bias_down_column_4_peak_vchannel","Avg BCAL peak vs channel;channel ID;peak",386,-1,1539);	

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

	bcal_peak_vevent->SetCanExtend(TH1::kXaxis);
	
	up_peak_vevent->SetCanExtend(TH1::kXaxis);
	down_peak_vevent->SetCanExtend(TH1::kXaxis);
	
	column1_peak_vevent->SetCanExtend(TH1::kXaxis);
	column2_peak_vevent->SetCanExtend(TH1::kXaxis);
	column3_peak_vevent->SetCanExtend(TH1::kXaxis);
	column4_peak_vevent->SetCanExtend(TH1::kXaxis);

	column1_up_peak_vevent->SetCanExtend(TH1::kXaxis);
	column2_up_peak_vevent->SetCanExtend(TH1::kXaxis);
	column3_up_peak_vevent->SetCanExtend(TH1::kXaxis);
	column4_up_peak_vevent->SetCanExtend(TH1::kXaxis);
	column1_down_peak_vevent->SetCanExtend(TH1::kXaxis);
	column2_down_peak_vevent->SetCanExtend(TH1::kXaxis);
	column3_down_peak_vevent->SetCanExtend(TH1::kXaxis);
	column4_down_peak_vevent->SetCanExtend(TH1::kXaxis);
	
	//////////////////////////////////////////////////////////////////////
#else
	bcal_peak_vevent->SetBit(TH1::kCanRebin);
	
	up_peak_vevent->SetBit(TH1::kCanRebin);
	down_peak_vevent->SetBit(TH1::kCanRebin);
	
	column1_peak_vevent->SetBit(TH1::kCanRebin);
	column2_peak_vevent->SetBit(TH1::kCanRebin);
	column3_peak_vevent->SetBit(TH1::kCanRebin);
	column4_peak_vevent->SetBit(TH1::kCanRebin);

	column1_up_peak_vevent->SetBit(TH1::kCanRebin);
	column2_up_peak_vevent->SetBit(TH1::kCanRebin);
	column3_up_peak_vevent->SetBit(TH1::kCanRebin);
	column4_up_peak_vevent->SetBit(TH1::kCanRebin);
	column1_down_peak_vevent->SetBit(TH1::kCanRebin);
	column2_down_peak_vevent->SetBit(TH1::kCanRebin);
	column3_down_peak_vevent->SetBit(TH1::kCanRebin);
	column4_down_peak_vevent->SetBit(TH1::kCanRebin);
	
	/////////////////////////////////////////////////////////
#endif

	h2_ledboth_Aall_vs_event = new TProfile("h2_ledboth_Aall_vs_event", "LED uboth - Aup and Adown vs event", 20000,0,200000000);
	h2_ledboth_sector_vs_event = new TProfile("h2_ledboth_sector_vs_event", "LED both - sector vs event", 20000,0,200000000);

	h1_ledup_sector = new TH1I("h1_ledup_sector", "LED up - sector", 50,0,5);
	h2_ledup_z_vs_cellid = new TH2I("h2_ledup_z_vs_cellid", "LED up - z vs Chan ID", 800,0,800,500,-100,400);
	h1_ledup_sector_config = new TH1I("h1_ledup_sector_config", "LED up - sector -config", 5,0,5);
	h1_ledup_z_all = new TH1I("h1_ledup_z_all", "LED up - z all channels", 500,-100,400);
	h1_ledup_Tdiff_all = new TH1I("h1_ledup_Tdiff_all", "LED up - Tdiff all channels", 400,-50,50);
	h1_ledup_Tup_all = new TH1I("h1_ledup_Tup_all", "LED up - Tup all channels", 410,0,410);
	h1_ledup_Tdown_all = new TH1I("h1_ledup_Tdown_all", "LED up - Tdown all channels", 410,0,410);
	h1_ledup_Aup_all = new TH1I("h1_ledup_Aup_all", "LED up - Aup all channels", 410,0,4100);
	h1_ledup_Adown_all = new TH1I("h1_ledup_Adown_all", "LED up - Adown all channels", 410,0,4100);
	h2_ledup_Aup_vs_z = new TH2I("h2_ledup_Aup_vs_z", "LED up - Aup vs z", 100,-100,400,410,0,4100);
	h2_ledup_Adown_vs_z = new TH2I("h2_ledup_Adown_vs_z", "LED up - Adown vs z", 100,-100,400,410,0,4100);
	h2_ledup_Aup_vs_event = new TProfile("h2_ledup_Aup_vs_event", "LED up - Aup vs event", 20000,0,200000000);
	h2_ledup_Adown_vs_event = new TProfile("h2_ledup_Adown_vs_event", "LED up - Adown vs event", 20000,0,200000000);
	h2_ledup_Aall_vs_event = new TProfile("h2_ledup_Aall_vs_event", "LED up - Aup and Adown vs event", 20000,0,200000000);
	h2_ledup_sector_vs_event = new TProfile("h2_ledup_sector_vs_event", "LED up - sector vs event", 20000,0,200000000);

	h1_leddown_sector = new TH1I("h1_leddown_sector", "LED down - sector", 50,0,5);
	h2_leddown_z_vs_cellid = new TH2I("h2_leddown_z_vs_cellid", "LED down - z vs Chan ID", 800,0,800,500,-100,400);
	h1_leddown_sector_config = new TH1I("h1_leddown_sector_config", "LED down - sector -config", 5,0,5);
	h1_leddown_z_all = new TH1I("h1_leddown_z_all", "LED down - z all channels", 500,-100,400);
	h1_leddown_Tdiff_all = new TH1I("h1_leddown_Tdiff_all", "LED down - Tdiff all channels", 400,-50,50);
	h1_leddown_Tup_all = new TH1I("h1_leddown_Tup_all", "LED down - Tup all channels", 410,0,410);
	h1_leddown_Tdown_all = new TH1I("h1_leddown_Tdown_all", "LED down - Tdown all channels", 410,0,410);
	h1_leddown_Aup_all = new TH1I("h1_leddown_Aup_all", "LED down - Aup all channels", 410,0,4100);
	h1_leddown_Adown_all = new TH1I("h1_leddown_Adown_all", "LED down - Adown all channels", 410,0,4100);
	h2_leddown_Aup_vs_z = new TH2I("h2_leddown_Aup_vs_z", "LED down - Aup vs z", 100,-100,400,410,0,4100);
	h2_leddown_Adown_vs_z = new TH2I("h2_leddown_Adown_vs_z", "LED down - Adown vs z", 100,-100,400,410,0,4100);
	h2_leddown_Aup_vs_event = new TProfile("h2_leddown_Aup_vs_event", "LED down - Aup vs event", 20000,0,200000000);
	h2_leddown_Adown_vs_event = new TProfile("h2_leddown_Adown_vs_event", "LED down - Adown vs event", 20000,0,200000000);
	h2_leddown_Aall_vs_event = new TProfile("h2_leddown_Aall_vs_event", "LED down - Aup and Adown vs event", 20000,0,200000000);
	h2_leddown_sector_vs_event = new TProfile("h2_leddown_sector_vs_event", "LED down - sector vs event", 20000,0,200000000);

	// back to main dir
	main->cd();
	
	// unlock
	japp->RootUnLock();
	
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED::evnt(JEventLoop *loop, uint64_t eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	
	int chcounter[1536] = { 0 } ;
	  
	
	vector<const DBCALDigiHit*> bcaldigihits;
	
	vector <const DBCALHit*> hitVector;

	vector<const DBCALHit*> dbcalhits;
	vector<const DBCALPoint*> dbcalpoints;
	
	
	bool LED_US=0, LED_DS=0;
	
	const DL1Trigger *trig = NULL;
	try {
		loop->GetSingle(trig);
	} catch (...) {}
	if (trig) {

		if (trig->trig_mask){
			// GTP tigger
			//GTPtrig++;
		}
		if (trig->fp_trig_mask){
			// Front panel trigger
			//FPtrig++;
		}
		if (trig->trig_mask && trig->fp_trig_mask){
			// Both GTP and front panel trigger
			//FPGTPtrig++;
		}
		if (trig->trig_mask & 0x1){
			// Cosmic trigger fired
			//trigCosmic++;
		}
		if (trig->fp_trig_mask & 0x100){
			// Upstream LED trigger fired
			//trigUS++;
			LED_US=1;
			
		}
		if (trig->fp_trig_mask & 0x200){
			// Downstream LED trigger fired
			//trigDS++;
			LED_DS=1;
		
		}
	} else {
		//NOtrig++;
	}	
	// Lock ROOT
	japp->RootWriteLock();

	float ledup_sector = 0;
	int ledup_sector_int = 0;
	float ledup_mean = 0;
	int ledup_events = 0;
	float leddown_sector = 0;
	int leddown_sector_int = 0;
	float leddown_mean = 0;
	int leddown_events = 0;
	
	if (LED_US || LED_DS) {
		

		loop->Get(dbcalhits);
		loop->Get(bcaldigihits);
		loop->Get(dbcalpoints);
		
	        // float apedsubtime[1536] = { 0. };
	        int apedsubpeak[1536] = { 0 };
		     
		for( unsigned int i=0; i<dbcalpoints.size(); i++) {

			int module = dbcalpoints[i]->module();
			int layer = dbcalpoints[i]->layer();
			int sector = dbcalpoints[i]->sector();
			int cell_id = (module-1)*16 + (layer-1)*4 + sector-1;
			float z =  dbcalpoints[i]->z();
			dbcalpoints[i]->Get(hitVector);
			const DBCALHit *Hit1 = hitVector[0];
			const DBCALHit *Hit2 = hitVector[1];
			// float Aup =  dbcalpoints[i]->E_US();
			// float Adown =  dbcalpoints[i]->E_DS();
			float Aup = 0;
			float Adown = 0;
			float Tup = 0;
			float Tdown = 0;
			float Tdiff = 0;
			if (Hit1->end == DBCALGeometry::kUpstream && Hit2->end == DBCALGeometry::kDownstream) {	  
			  Aup = Hit1->pulse_peak;
			  Adown = Hit2->pulse_peak;	  
			  Tup = Hit1->t_raw;
			  Tdown = Hit2->t_raw;
			  Tdiff = Tdown - Tup;
			}
			else if (Hit2->end == DBCALGeometry::kUpstream && Hit1->end == DBCALGeometry::kDownstream){
			  Aup = Hit2->pulse_peak;
			  Adown = Hit1->pulse_peak;
			  Tup = Hit2->t_raw;
			  Tdown = Hit1->t_raw;
			  Tdiff = Tdown - Tup;
			}

			// fill histograms for all channels
			if (LED_US) {
			  h1_ledup_sector->Fill(sector);
			  h1_ledup_z_all->Fill(z);;
			  h2_ledup_z_vs_cellid->Fill(cell_id,z);
			  h1_ledup_Tdiff_all->Fill(Tdiff);;
			  h1_ledup_Tup_all->Fill(Tup);;
			  h1_ledup_Tdown_all->Fill(Tdown);;
			  h1_ledup_Aup_all->Fill(Aup);;
			  h1_ledup_Adown_all->Fill(Adown);;
			  h2_ledup_Aup_vs_z->Fill(z,Aup);
			  h2_ledup_Adown_vs_z->Fill(z,Adown);
			  h2_ledup_Aup_vs_event->Fill(eventnumber,Aup);
			  h2_ledup_Adown_vs_event->Fill(eventnumber,Adown);
			  // cout << "US - eventnumber=" << eventnumber << " Aup=" << Aup << " Aup=" << Aup << " Adown=" << Adown << " Adown=" << Adown << endl;
			}
			else if (LED_DS) {
			  h1_leddown_sector->Fill(sector);
			  h1_leddown_z_all->Fill(z);
			  h2_leddown_z_vs_cellid->Fill(cell_id,z);
			  h1_leddown_Tdiff_all->Fill(Tdiff);;
			  h1_leddown_Tup_all->Fill(Tup);;
			  h1_leddown_Tdown_all->Fill(Tdown);;
			  h1_leddown_Aup_all->Fill(Aup);;
			  h1_leddown_Adown_all->Fill(Adown);;
			  h2_leddown_Aup_vs_z->Fill(z,Aup);
			  h2_leddown_Adown_vs_z->Fill(z,Adown);
			  h2_leddown_Aup_vs_event->Fill(eventnumber,Aup);
			  h2_leddown_Adown_vs_event->Fill(eventnumber,Adown);
			  // cout << "DS - eventnumber=" << eventnumber << " Aup=" << Aup << " Aup=" << Aup << " Adown=" << Adown << " Adown=" << Adown << endl;
			}

			// make cuts on z for all hits
			// if ( LED_DS && (z>280 && z<400)) {
			// cut on Tdiff instead
			if ( LED_DS && (Tdiff>-30 && Tdiff<-15)) {
				apedsubpeak[cell_id] = Adown;
				chcounter[cell_id]++;
				apedsubpeak[cell_id+768] = Aup;
				chcounter[cell_id+768]++;
				h2_ledboth_Aall_vs_event->Fill(eventnumber,Aup);
				h2_ledboth_Aall_vs_event->Fill(eventnumber,Adown);
				h2_leddown_Aall_vs_event->Fill(eventnumber,Aup);
				h2_leddown_Aall_vs_event->Fill(eventnumber,Adown);
				h2_leddown_sector_vs_event->Fill(eventnumber,sector);
				h2_ledboth_sector_vs_event->Fill(eventnumber,sector);
				
				// compute sums
				leddown_sector += sector;
				leddown_mean += Aup + Adown;
				leddown_events++;
				
			} // if condition on Tdiff

			// if (LED_US && (z>-60 && z<0)) {
			if (LED_US && (Tdiff>15 && Tdiff<30)) {
				apedsubpeak[cell_id] = Adown;
				chcounter[cell_id]++;
				apedsubpeak[cell_id+768] = Aup;
				chcounter[cell_id+768]++;
				h2_ledboth_Aall_vs_event->Fill(eventnumber,Aup);
				h2_ledboth_Aall_vs_event->Fill(eventnumber,Adown);
				h2_ledup_Aall_vs_event->Fill(eventnumber,Aup);
				h2_ledup_Aall_vs_event->Fill(eventnumber,Adown);
				h2_ledup_sector_vs_event->Fill(eventnumber,sector);
				h2_ledboth_sector_vs_event->Fill(eventnumber,sector);
				
				// compute sums
				ledup_sector += sector;
				ledup_mean += Aup + Adown;
				ledup_events++;
			} // if condition on z

		}//loop over bcalhits

	 //  compute averages
	 ledup_sector_int = ledup_events > 0? ledup_sector/ledup_events + 0.5 : 0;
	 ledup_mean = ledup_events > 0? ledup_mean/(2*ledup_events): 0;
	 leddown_sector_int = leddown_events > 0? leddown_sector/leddown_events + 0.5: 0;
	 leddown_mean = leddown_events > 0? leddown_mean/(2*leddown_events) : 0;
	 if (LED_US) h1_ledup_sector_config->Fill(ledup_sector_int);
	 if (LED_DS) h1_leddown_sector_config->Fill(leddown_sector_int);

	 // cout << " ledup_evens=" << ledup_events << " ledup_sector=" << ledup_sector << " ledup_mean=" << ledup_mean << endl;
	 // cout << " leddown_evens=" << leddown_events << " ledown_sector=" << leddown_sector << " leddown_mean=" << leddown_mean << endl << endl;

			
	 // float sector_delta=0.2;
	for (int chid = 0; chid < 1536; chid++)  {
	  // if (chcounter[chid] > 1) continue;
	      if (chcounter[chid] != 1) continue;
	      bcal_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				if (LED_US) {
				    up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				    if (ledup_sector_int == 1 && chid%4+1 == 1) {column1_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (ledup_sector_int == 2 && chid%4+1 == 2) {column2_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (ledup_sector_int == 3 && chid%4+1 == 3) {column3_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (ledup_sector_int == 4 && chid%4+1 == 4) {column4_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				   }
				    
				else if (LED_DS) {
				    down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				    if (leddown_sector_int == 1 && chid%4+1 == 1) {column1_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (leddown_sector_int == 2 && chid%4+1 == 2) {column2_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (leddown_sector_int == 3 && chid%4+1 == 3) {column3_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (leddown_sector_int == 4 && chid%4+1 == 4) {column4_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    }
				    
				if (ledup_sector_int == 1 && chid%4+1 == 1) {column1_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (ledup_sector_int == 2 && chid%4+1 == 2) {column2_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (ledup_sector_int == 3 && chid%4+1 == 3) {column3_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (ledup_sector_int == 4 && chid%4+1 == 4) {column4_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
	    }//loop over bcalhits

		   //Deduce LED pulsing configuration based on average pulse peak in BCAL, each side & each column then fill correponding profile.
		 double column1up = 0;
 		 double column2up = 0;
		 double column3up = 0;
		 double column4up = 0;

		 double column1down = 0;
		 double column2down = 0;
 		 double column3down = 0;
		 double column4down = 0;

		 column1up = column1_up_peak_vevent->GetBinContent(column1_up_peak_vevent->FindBin(eventnumber));
 		 column2up = column2_up_peak_vevent->GetBinContent(column2_up_peak_vevent->FindBin(eventnumber));
		 column3up = column3_up_peak_vevent->GetBinContent(column3_up_peak_vevent->FindBin(eventnumber));
		 column4up = column4_up_peak_vevent->GetBinContent(column4_up_peak_vevent->FindBin(eventnumber));

		 column1down = column1_down_peak_vevent->GetBinContent(column1_down_peak_vevent->FindBin(eventnumber));
		 column2down = column2_down_peak_vevent->GetBinContent(column2_down_peak_vevent->FindBin(eventnumber));
 		 column3down = column3_down_peak_vevent->GetBinContent(column3_down_peak_vevent->FindBin(eventnumber));
		 column4down = column4_down_peak_vevent->GetBinContent(column4_down_peak_vevent->FindBin(eventnumber));

		 // cout << " adccount1=" << adccount1 << " column1up=" << column1up << " column2up=" << column2up << " column3up=" << column3up << " column4up=" << column4up << endl;
		 // cout << " adccount1=" << adccount1 << " column1down=" << column1down << " column2down=" << column2down << " column3down=" << column3down << " column4down=" << column4down << endl << endl;

		 // Now categorize according to configuration

	for (int chid = 0; chid < 1536; chid++)  {
	  // if (chcounter[chid] > 1) continue;
	      if (chcounter[chid] != 1) continue;
				if (LED_US) {
				    if (ledup_sector_int == 1 && chid%4+1 == 1) {
				      if (column1up < adccount1) {
					low_up_1->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column1_down_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column1_up_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_up_1->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column1_down_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column1_up_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (ledup_sector_int == 2 && chid%4+1 == 2) {
				      if (column2up < adccount1) {
					low_up_2->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column2_down_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column2_up_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_up_2->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column2_down_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column2_up_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (ledup_sector_int == 3 && chid%4+1 == 3)  {
				      if (column3up < adccount1) {
					low_up_3->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column3_down_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column3_up_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_up_3->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column3_down_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column3_up_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (ledup_sector_int == 4 && chid%4+1 == 4) {
				      if (column4up < adccount1) {
					low_up_4->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column4_down_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column4_up_peak_vevent2->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_up_4->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column4_down_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column4_up_peak_vevent4->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				}	    
				else if (LED_DS) {
				    if (leddown_sector_int == 1 && chid%4+1 == 1) {
				      if (column1down < adccount1) {
					low_down_1->Fill(chid,apedsubpeak[chid]);
					// if (eventnumber >31500000 && eventnumber <46000000) cout << " chid=" << chid << " column1down=" << column1down << " apedsub=" << apedsubpeak[chid] << endl;
					if (chid < 768) {
					  column1_down_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column1_up_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_down_1->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column1_down_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column1_up_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (leddown_sector_int == 2 && chid%4+1 == 2) {
				      if (column2down < adccount1) {
					low_down_2->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column2_down_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column2_up_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_down_2->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column2_down_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column2_up_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (leddown_sector_int == 3 && chid%4+1 == 3)  {
				      if (column3down < adccount1) {
					low_down_3->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column3_down_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column3_up_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_down_3->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column3_down_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column3_up_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				    else if (leddown_sector_int == 4 && chid%4+1 == 4) {
				      if (column4down < adccount1) {
					low_down_4->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column4_down_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column4_up_peak_vevent1->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				      else {
					high_down_4->Fill(chid,apedsubpeak[chid]);
					if (chid < 768) {
					  column4_down_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					  }
					 else if  (chid > 767) {
					   column4_up_peak_vevent3->Fill(eventnumber,apedsubpeak[chid]);
					 }
				      }
				    }
				}
			
	}  //loop over bcalhits

	}//if LEDUP || LEDDOWN    
	// Unlock ROOT
	japp->RootUnLock();
	

    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED::erun(void) {
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.


/*	printf("\nTrigger statistics");
	printf("------------------------\n");
	printf("%20s: %10i\n","no triggers",NOtrig);
	printf("%20s: %10i\n","Front Panel",FPtrig);
	printf("%20s: %10i\n","GTP",GTPtrig);
	printf("%20s: %10i\n","FP && GTP",FPGTPtrig);
	printf("%20s: %10i\n","US LED",trigUS);
	printf("%20s: %10i\n","DS LED",trigDS);
	printf("%20s: %10i\n","BCAL",trigCosmic);
	ledcounter = low_down_1_counter + low_down_2_counter + low_down_3_counter + low_down_4_counter + low_up_1_counter + low_up_2_counter + low_up_3_counter + low_up_4_counter + high_down_1_counter + high_down_2_counter + high_down_3_counter + high_down_4_counter + high_up_1_counter + high_up_2_counter + high_up_3_counter + high_up_4_counter + unidentified;
	printf("%20s: %10i\n","low_down_1_counter",low_down_1_counter);
	printf("%20s: %10i\n","low_down_2_counter",low_down_2_counter);
	printf("%20s: %10i\n","low_down_3_counter",low_down_3_counter);
	printf("%20s: %10i\n","low_down_4_counter",low_down_4_counter);

	printf("%20s: %10i\n","low_up_1_counter",low_up_1_counter);
	printf("%20s: %10i\n","low_up_2_counter",low_up_2_counter);
	printf("%20s: %10i\n","low_up_3_counter",low_up_3_counter);
	printf("%20s: %10i\n","low_up_4_counter",low_up_4_counter);

	printf("%20s: %10i\n","high_down_1_counter",high_down_1_counter);
	printf("%20s: %10i\n","high_down_2_counter",high_down_2_counter);
	printf("%20s: %10i\n","high_down_3_counter",high_down_3_counter);
	printf("%20s: %10i\n","high_down_4_counter",high_down_4_counter);

	printf("%20s: %10i\n","high_up_1_counter",high_up_1_counter);
	printf("%20s: %10i\n","high_up_2_counter",high_up_2_counter);
	printf("%20s: %10i\n","high_up_3_counter",high_up_3_counter);
	printf("%20s: %10i\n","high_up_4_counter",high_up_4_counter);
	
	printf("%20s: %10i\n","Unidentified",unidentified);*/
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED::fini(void) {
	// Called before program exit after event processing is finished.

// 	//Write mean pulse peak to output file
	ofstream foutlowdown ; foutlowdown.open("LED_lowbias_downstream.txt");
	ofstream foutlowup; foutlowup.open("LED_lowbias_upstream.txt");
	ofstream fouthighdown; fouthighdown.open("LED_highbias_downstream.txt");
	ofstream fouthighup; fouthighup.open("LED_highbias_upstream.txt");

	for(int k=0 ;k < 768;k += 4)
	{
	double lowdownmean1down = low_down_1->GetBinContent(low_down_1->FindBin(k));
	double lowdownmean1up = low_down_1->GetBinContent(low_down_1->FindBin(768+k));
	double lowdownmean2down = low_down_2->GetBinContent(low_down_2->FindBin(k+1));
	double lowdownmean2up = low_down_2->GetBinContent(low_down_2->FindBin(769+k));
	double lowdownmean3down = low_down_3->GetBinContent(low_down_3->FindBin(k+2));
	double lowdownmean3up = low_down_3->GetBinContent(low_down_3->FindBin(770+k));
	double lowdownmean4down = low_down_4->GetBinContent(low_down_4->FindBin(k+3));
	double lowdownmean4up = low_down_4->GetBinContent(low_down_4->FindBin(771+k));
	
	double lowupmean1down = low_up_1->GetBinContent(low_up_1->FindBin(k));
	double lowupmean1up = low_up_1->GetBinContent(low_up_1->FindBin(768+k));
	double lowupmean2down = low_up_2->GetBinContent(low_up_2->FindBin(k+1));
	double lowupmean2up = low_up_2->GetBinContent(low_up_2->FindBin(769+k));
	double lowupmean3down = low_up_3->GetBinContent(low_up_3->FindBin(k+2));
	double lowupmean3up = low_up_3->GetBinContent(low_up_3->FindBin(770+k));
	double lowupmean4down = low_up_4->GetBinContent(low_up_4->FindBin(k+3));
	double lowupmean4up = low_up_4->GetBinContent(low_up_4->FindBin(771+k));

	double highdownmean1down = high_down_1->GetBinContent(high_down_1->FindBin(k));
	double highdownmean1up = high_down_1->GetBinContent(high_down_1->FindBin(768+k));
	double highdownmean2down = high_down_2->GetBinContent(high_down_2->FindBin(k+1));
	double highdownmean2up = high_down_2->GetBinContent(high_down_2->FindBin(769+k));
	double highdownmean3down = high_down_3->GetBinContent(high_down_3->FindBin(k+2));
	double highdownmean3up = high_down_3->GetBinContent(high_down_3->FindBin(770+k));
	double highdownmean4down = high_down_4->GetBinContent(high_down_4->FindBin(k+3));
	double highdownmean4up = high_down_4->GetBinContent(high_down_4->FindBin(771+k));
	
	double highupmean1down = high_up_1->GetBinContent(high_up_1->FindBin(k));
	double highupmean1up = high_up_1->GetBinContent(high_up_1->FindBin(768+k));
	double highupmean2down = high_up_2->GetBinContent(high_up_2->FindBin(k+1));
	double highupmean2up = high_up_2->GetBinContent(high_up_2->FindBin(769+k));
	double highupmean3down = high_up_3->GetBinContent(high_up_3->FindBin(k+2));
	double highupmean3up = high_up_3->GetBinContent(high_up_3->FindBin(770+k));
	double highupmean4down = high_up_4->GetBinContent(high_up_4->FindBin(k+3));
	double highupmean4up = high_up_4->GetBinContent(high_up_4->FindBin(771+k));
	
	//TString sep = "        ";
	foutlowdown << lowdownmean1down << endl << lowdownmean1up << endl << lowdownmean2down << endl << lowdownmean2up << endl << lowdownmean3down << endl << lowdownmean3up << endl << lowdownmean4down << endl << lowdownmean4up << endl;
	foutlowup << lowupmean1down << endl << lowupmean1up << endl << lowupmean2down << endl << lowupmean2up << endl << lowupmean3down << endl << lowupmean3up << endl << lowupmean4down << endl << lowupmean4up << endl;
	fouthighdown << highdownmean1down << endl << highdownmean1up << endl << highdownmean2down << endl << highdownmean2up << endl << highdownmean3down << endl << highdownmean3up << endl << highdownmean4down << endl << highdownmean4up << endl;
	fouthighup << highupmean1down << endl << highupmean1up << endl << highupmean2down << endl << highupmean2up << endl << highupmean3down << endl << highupmean3up << endl << highupmean4down << endl << highupmean4up << endl;

	}

	foutlowdown.close();
	foutlowup.close();
	fouthighdown.close();
	fouthighup.close();


return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
