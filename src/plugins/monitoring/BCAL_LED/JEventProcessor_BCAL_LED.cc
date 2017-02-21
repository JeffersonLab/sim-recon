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


	maxnumberofevents=1000;
	nbins=1002;//Extendable histograms with a single event per bin
	
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

	// back to main dir
	main->cd();
	
	// unlock
	japp->RootUnLock();
	
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	
	if(runnumber > 9999 && runnumber < 20000)//Spring 2016 run period
	{
	adccount1 = 1100;
	adccount2 = 1400;
	adccount3 = 1800;
	}
	else if (runnumber > 29999 && runnumber < 40000)//Spring 2017 run period
	{
        adccount1=600;
        adccount2=900;
        adccount3=1400;
	}
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
	

	vector<const DBCALHit*> dbcalhits;
	
	
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
	
	if (LED_US || LED_DS) {
		

		loop->Get(dbcalhits);
		loop->Get(bcaldigihits);
		
		     int apedsubpeak[1536] = { 0 };
		     
		     int cellsector[1536] =  { 0 };
		     int cellend[1536] =  { 0 };

		     
		for( unsigned int i=0; i<dbcalhits.size(); i++) {
                        const DBCALHit *bcalhit = dbcalhits[i];
			const DBCALDigiHit *bcaldigihit = NULL;
			bcalhit->GetSingle(bcaldigihit);

			int module = bcalhit->module;
			int layer = bcalhit->layer;
			int sector = bcalhit->sector;
		        int end = bcalhit->end;
			int cell_id = -1;
				
				
				if(bcalhit->end == DBCALGeometry::kDownstream)
				{
				cell_id = (module-1)*16 + (layer-1)*4 + sector-1; //has a range of 768 channels enough for one side.

				}//if downstream cell id
				if(bcalhit->end == DBCALGeometry::kUpstream)
				{
				cell_id = 768 + (module-1)*16 + (layer-1)*4 + sector-1;
				}
				
				apedsubpeak[cell_id] = bcaldigihit->pulse_peak - (int) bcaldigihit->pedestal / bcaldigihit->nsamples_pedestal;
				chcounter[cell_id]++;
				
				cellsector[cell_id] = sector;
				cellend[cell_id] = end;

        }//loop over bcalhits
				
		// Lock ROOT
		japp->RootWriteLock();

	for (int chid = 0; chid < 1536; chid++)
	    {
	      if (chcounter[chid] > 1) continue;
	      bcal_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				if (cellend[chid] == DBCALGeometry::kUpstream)
				   {up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				    if (chid%4 == 0) {column1_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 1) {column2_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 2) {column3_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 3) {column4_up_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				   }
				    
				else if (cellend[chid] == DBCALGeometry::kDownstream) 
				    {down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);
				    if (chid%4 == 0) {column1_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 1) {column2_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 2) {column3_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    else if (chid%4 == 3) {column4_down_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				    }
				    
				if (cellsector[chid] == 1) {column1_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (cellsector[chid] == 2) {column2_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (cellsector[chid] == 3) {column3_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
				else if (cellsector[chid] == 4) {column4_peak_vevent->Fill(eventnumber,apedsubpeak[chid]);}
		}//loop over bcalhits

		   //Deduce LED pulsing configuration based on average pulse peak in BCAL, each side & each column then fill correponding profile.
		 double column1up = column1_up_peak_vevent->GetBinContent(column1_up_peak_vevent->FindBin(eventnumber));
 		 double column2up = column2_up_peak_vevent->GetBinContent(column2_up_peak_vevent->FindBin(eventnumber));
		 double column3up = column3_up_peak_vevent->GetBinContent(column3_up_peak_vevent->FindBin(eventnumber));
		 double column4up = column4_up_peak_vevent->GetBinContent(column4_up_peak_vevent->FindBin(eventnumber));

		 double column1down = column1_down_peak_vevent->GetBinContent(column1_down_peak_vevent->FindBin(eventnumber));
		 double column2down = column2_down_peak_vevent->GetBinContent(column2_down_peak_vevent->FindBin(eventnumber));
 		 double column3down = column3_down_peak_vevent->GetBinContent(column3_down_peak_vevent->FindBin(eventnumber));
		 double column4down = column4_down_peak_vevent->GetBinContent(column4_down_peak_vevent->FindBin(eventnumber));
		 
			if      (adccount1 < column1up && column1up < adccount2 && column1up > column1down && column1up > column2up && column1down > column2down && column1up > column3up && column1down > column3down && column1up > column4up && column1down > column4down)
			    {
			      //column = 1;
			    for(int k=0 ;k < 1536;k++) 
			    {
			      if (chcounter[k] > 1) continue;
			      if (k%4 == 0 && apedsubpeak[k] > 0) {low_down_1->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column1_down_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column1_up_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_down_1_counter++;
			    }
			else if (adccount1 < column2up && column2up < adccount2 && column2up > column2down && column2up > column1up && column2down > column1down && column2up > column3up && column2down > column3down && column2up > column4up && column2down > column4down)
			    {//column = 2;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 1 && apedsubpeak[k] > 0) {low_down_2->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column2_down_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column2_up_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_down_2_counter++;
			    }
			else if (adccount1 < column3up && column3up < adccount2 && column3up > column3down && column3up > column1up && column3down > column1down && column3up > column2up && column3down > column2down && column3up > column4up && column3down > column4down)
			    {//column = 3;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 2 && apedsubpeak[k] > 0) {low_down_3->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column3_down_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column3_up_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_down_3_counter++;
			    }
			else if (adccount1 < column4up && column4up < adccount2 && column4up > column4down && column4up > column1up && column4down > column1down && column4up > column2up && column4down > column2down && column4up > column3up && column4down > column3down)
			    {//column = 4;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 3 && apedsubpeak[k] > 0) {low_down_4->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column4_down_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column4_up_peak_vevent1->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_down_4_counter++;
			    }
			   
		    
			else if (adccount1 < column1down && column1down < adccount2 && column1down > column1up && column1up > column2up && column1down > column2down && column1up > column3up && column1down > column3down && column1up > column4up && column1down > column4down)
			    {//column = 1;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 0 && apedsubpeak[k] > 0) {low_up_1->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column1_down_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column1_up_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_up_1_counter++;			    
			    }
			    
			else if (adccount1 < column2down && column2down < adccount2 && column2down > column2up && column2up > column1up && column2down > column1down && column2up > column3up && column2down > column3down && column2up > column4up && column2down > column4down)
			    {//column = 2;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 1 && apedsubpeak[k] > 0) {low_up_2->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column2_down_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column2_up_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_up_2_counter++;
			    }
			else if (adccount1 < column3down && column3down < adccount2 && column3down > column3up && column3up > column1up && column3down > column1down && column3up > column2up && column3down > column2down && column3up > column4up && column3down > column4down)
			    {//column = 3;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 2 && apedsubpeak[k] > 0) {low_up_3->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column3_down_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column3_up_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_up_3_counter++;
			    }
			else if (adccount1 < column4down && column4down < adccount2 && column4down > column4up && column4up > column1up && column4down > column1down && column4up > column2up && column4down > column2down && column4up > column3up && column4down > column3down)
			    {//column = 4;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 3 && apedsubpeak[k] > 0) {low_up_4->Fill(k, apedsubpeak[k]);
					    if (k < 768) {column4_down_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column4_up_peak_vevent2->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //low_up_4_counter++;
			    }
			 
		      
		
			else if (column1up > adccount3 && column1up > column1down && column1up > column2up && column1down > column2down && column1up > column3up && column1down > column3down && column1up > column4up && column1down > column4down)
			    {//column = 1;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 0 && apedsubpeak[k] > 0) {high_down_1->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column1_down_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column1_up_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //high_down_1_counter++;
			    }
			else if (column2up > adccount3 && column2up > column2down && column2up > column1up && column2down > column1down && column2up > column3up && column2down > column3down && column2up > column4up && column2down > column4down)
			    {//column = 2;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 1 && apedsubpeak[k] > 0) {high_down_2->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column2_down_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column2_up_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //high_down_2_counter++;
			    }
			else if (column3up > adccount3 && column3up > column3down && column3up > column1up && column3down > column1down && column3up > column2up && column3down > column2down && column3up > column4up && column3down > column4down)
			    {//column = 3;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 2 && apedsubpeak[k] > 0) {high_down_3->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column3_down_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column3_up_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //high_down_3_counter++;
			    }
			else if (column4up > adccount3 && column4up > column4down && column4up > column1up && column4down > column1down && column4up > column2up && column4down > column2down && column4up > column3up && column4down > column3down)
			    {//column = 4;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 3 && apedsubpeak[k] > 0) {high_down_4->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column4_down_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column4_up_peak_vevent3->Fill(eventnumber,apedsubpeak[k]);
								}
					    }
			    }
			    //high_down_4_counter++;
			    }
			  

			else if (column1down > adccount3 && column1down > column1up && column1up > column2up && column1down > column2down && column1up > column3up && column1down > column3down && column1up > column4up && column1down > column4down)
			    {//column = 1;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 0 && apedsubpeak[k] > 0) {high_up_1->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column1_down_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column1_up_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
								}
					   }
			    }
			    //high_up_1_counter++;
			    }
			else if (column2down > adccount3 && column2down > column2up && column2up > column1up && column2down > column1down && column2up > column3up && column2down > column3down && column2up > column4up && column2down > column4down)
			    {//column = 2;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 1 && apedsubpeak[k] > 0) {high_up_2->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column2_down_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column2_up_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
								}
					   }
			    }
			    //high_up_2_counter++;
			    }
			else if (column3down > adccount3 && column3down > column3up && column3up > column1up && column3down > column1down && column3up > column2up && column3down > column2down && column3up > column4up && column3down > column4down)
			    {//column = 3;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 2 && apedsubpeak[k] > 0) {high_up_3->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column3_down_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column3_up_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
								}
					   }
			    }
			    //high_up_3_counter++;
			    }
			else if (column4down > adccount3 && column4down > column4up && column4up > column1up && column4down > column1down && column4up > column2up && column4down > column2down && column4up > column3up && column4down > column3down)
			    {//column = 4;
			    for(int k=0 ;k < 1536;k++) 
			    {if (k%4 == 3 && apedsubpeak[k] > 0) {high_up_4->Fill(k, apedsubpeak[k]);
    					    if (k < 768) {column4_down_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
							  }
					    else if  (k > 767) {column4_up_peak_vevent4->Fill(eventnumber,apedsubpeak[k]);
								}
					   }
			    }
			    //high_up_4_counter++;
			    }
			//else {unidentified++;}    
		     // Unlock ROOT
		japp->RootUnLock();
	}//if LEDUP || LEDDOWN
	

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
