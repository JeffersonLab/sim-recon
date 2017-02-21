// $Id$
//
//    File: JEventProcessor_BCAL_LED_time.cc
//

#include <stdint.h>
#include <vector>
#include "TTree.h"
#include "JEventProcessor_BCAL_LED_time.h"
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


static TH1I* h1_ledup_z_all = NULL;
static TH2I* h2_ledup_z_vs_cellid = NULL;
static TH2I* h2_ledup_Eup_vs_cellid = NULL;
static TH2I* h2_ledup_Edown_vs_cellid = NULL;
static TH2I* h2_ledup_Eup_vs_z = NULL;
static TH2I* h2_ledup_Edown_vs_z = NULL;
static TH2I* h2_ledup_z_vs_event = NULL;
static TH2I* h2_ledup_Eup_vs_Edown = NULL;

static TH1I* h1_leddown_z_all = NULL;
static TH2I* h2_leddown_z_vs_cellid = NULL;
static TH2I* h2_leddown_Eup_vs_cellid = NULL;
static TH2I* h2_leddown_Edown_vs_cellid = NULL;
static TH2I* h2_leddown_Eup_vs_z = NULL;
static TH2I* h2_leddown_Edown_vs_z = NULL;
static TH2I* h2_leddown_z_vs_event = NULL;
static TH2I* h2_leddown_Eup_vs_Edown = NULL;

//all channels
static TProfile *bcal_peak_vevent = NULL;


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
		app->AddProcessor(new JEventProcessor_BCAL_LED_time());
	}

}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LED_time::JEventProcessor_BCAL_LED_time() {
}


//----------------------------------------------------------------------------------


JEventProcessor_BCAL_LED_time::~JEventProcessor_BCAL_LED_time() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_BCAL_LED_time::init(void) {
	
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
	adccount2 = 1400;
	adccount3 = 1800;
	
	maxnumberofevents=700000000.0;//Assuming 1Hz LED trigger, 300M for a beam run with 30KHz trigger and 700M for 70KHz
	//maxnumberofevents=10000.0;//using LED event conter
	nbins=24002;//Assuming 1Hz LED trigger, 10K for a beam run with 30KHz trigger and 24K for 70KHz
	//nbins=375002;//Based on cosmic run with 800Hz trigger and 1Hz LED trigger
	//nbins=3750002;//Based on cosmic run with 800Hz trigger and 10Hz LED trigger

	
	//overflow=0; underflow=0; negatives=0; zeros=0;

	// create root folder for bcal and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("bcalLED")->cd();
	gStyle->SetOptStat(111110);

	bcal_peak_vevent = new TProfile("bcal_peak_vevent","Avg BCAL peak vs event;event num;peak (all chan avg)",nbins,0.0,maxnumberofevents);

	h1_ledup_z_all = new TH1I("h1_ledup_z_all", "LED up - z all channels", 500,-100,400);
	h2_ledup_z_vs_cellid = new TH2I("h2_ledup_z_vs_cellid", "LED up - z vs Chan ID", 800,0,800,500,-100,400);
	h2_ledup_Eup_vs_cellid = new TH2I("h2_ledup_Eup_vs_cellid", "LED up - Eup vs Chan ID", 800,0,800,400,0,4);
	h2_ledup_Edown_vs_cellid = new TH2I("h2_ledup_Edown_vs_cellid", "LED up - Edown vs Chan ID", 800,0,800,400,0,4);
	h2_ledup_Eup_vs_z = new TH2I("h2_ledup_Eup_vs_z", "LED up - Eup vs z", 100,-100,400,400,0,4);
	h2_ledup_Edown_vs_z = new TH2I("h2_ledup_Edown_vs_z", "LED up - Edown vs z", 100,-100,400,400,0,4);
	h2_ledup_z_vs_event = new TH2I("h2_ledup_z_vs_event", "LED up - z vs event number", 1000,0,100000000,500,-100,400);
	h2_ledup_Eup_vs_Edown = new TH2I("h2_ledup_Eup_vs_Edown", "LED up - Eup vs Edown", 400,0,4,400,0,4);

	h1_leddown_z_all = new TH1I("h1_leddown_z_all", "LED down - z all channels", 500,-100,400);
	h2_leddown_z_vs_cellid = new TH2I("h2_leddown_z_vs_cellid", "LED down - z vs Chan ID", 800,0,800,500,-100,400);
	h2_leddown_Eup_vs_cellid = new TH2I("h2_leddown_Eup_vs_cellid", "LED down - Eup vs Chan ID", 800,0,800,400,0,4);
	h2_leddown_Edown_vs_cellid = new TH2I("h2_leddown_Edown_vs_cellid", "LED down - Edown vs Chan ID", 800,0,800,400,0,4);
	h2_leddown_Eup_vs_z = new TH2I("h2_leddown_Eup_vs_z", "LED down - Eup vs z", 100,-100,400,400,0,4);
	h2_leddown_Edown_vs_z = new TH2I("h2_leddown_Edown_vs_z", "LED down - Edown vs z", 100,-100,400,400,0,4);
	h2_leddown_z_vs_event = new TH2I("h2_leddown_z_vs_event", "LED down - z vs event number", 1000,0,100000000,500,-100,400);
	h2_leddown_Eup_vs_Edown = new TH2I("h2_leddown_Eup_vs_Edown", "LED up - Eup vs Edown", 400,0,4,400,0,4);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)

	bcal_peak_vevent->SetCanExtend(TH1::kXaxis);
	
	//////////////////////////////////////////////////////////////////////
#else
	bcal_peak_vevent->SetBit(TH1::kCanRebin);
	
	/////////////////////////////////////////////////////////
#endif

	// back to main dir
	main->cd();
	
	// unlock
	japp->RootUnLock();
	
	
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED_time::brun(JEventLoop *eventLoop, int32_t runnumber) {
	// This is called whenever the run number changes
	return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED_time::evnt(JEventLoop *loop, uint64_t eventnumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop-Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	  
	
	vector<const DBCALDigiHit*> bcaldigihits;
	
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
			cout << "Upstream LED_US=" << LED_US << endl;
			
		}
		if (trig->fp_trig_mask & 0x200){
			// Downstream LED trigger fired
			//trigDS++;
			LED_DS=1;;
			cout << "Downstream LED_DS=" << LED_DS << endl;
		
		}
	} else {
		//NOtrig++;
	}
	
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
	if (LED_US) {
	        local_eventnum++;
		
		loop->Get(dbcalhits);
		loop->Get(bcaldigihits);
		loop->Get(dbcalpoints);
		

		for( unsigned int i=0; i<dbcalpoints.size(); i++) {
		  int module = dbcalpoints[i]->module();
		  int layer = dbcalpoints[i]->layer();
		  int sector = dbcalpoints[i]->sector();
		  int cell_id = (module-1)*16 + (layer-1)*4 + sector-1;
		  float z =  dbcalpoints[i]->z();
		  float Eup =  dbcalpoints[i]->E_US();
		  float Edown =  dbcalpoints[i]->E_DS();
		  // cout << " Up i=" <<  i << " module=" << module << " layer=" << layer << " sector=" << sector << " z=" << z << endl;
		  h1_ledup_z_all->Fill(z);
		  h2_ledup_z_vs_cellid->Fill(cell_id,z);
		  h2_ledup_Eup_vs_cellid->Fill(cell_id,Eup);
		  h2_ledup_Edown_vs_cellid->Fill(cell_id,Edown);
		  h2_ledup_Eup_vs_z->Fill(z,Eup);
		  h2_ledup_Edown_vs_z->Fill(z,Edown);
		  // cout << " UP i=" << i << " eventnumber=" << eventnumber << " local_eventnum=" << local_eventnum << endl;
		  // h2_ledup_z_vs_event->Fill(local_eventnum,z);    
		  h2_ledup_z_vs_event->Fill(eventnumber,z);  
		  h2_ledup_Eup_vs_Edown->Fill(Edown,Eup);
		}
	}
	if (LED_DS) {
	        local_eventnum++;
		
		loop->Get(dbcalhits);
		loop->Get(bcaldigihits);
		loop->Get(dbcalpoints);
		

		for( unsigned int i=0; i<dbcalpoints.size(); i++) {
		  int module = dbcalpoints[i]->module();
		  int layer = dbcalpoints[i]->layer();
		  int sector = dbcalpoints[i]->sector();
		  int cell_id = (module-1)*16 + (layer-1)*4 + sector-1;
		  float z =  dbcalpoints[i]->z();
		  float Eup =  dbcalpoints[i]->E_US();
		  float Edown =  dbcalpoints[i]->E_DS();
		  // cout << " Down i=" <<  i << " module=" << module << " layer=" << layer << " sector=" << sector << " z=" << z << endl;
		  h1_leddown_z_all->Fill(z);
		  h2_leddown_z_vs_cellid->Fill(cell_id,z);
		  h2_leddown_Eup_vs_cellid->Fill(cell_id,Eup);
		  h2_leddown_Edown_vs_cellid->Fill(cell_id,Edown);
		  h2_leddown_Eup_vs_z->Fill(z,Eup);
		  h2_leddown_Edown_vs_z->Fill(z,Edown);
		  // cout << " Down i=" << i << " eventnumber=" << eventnumber << " local_eventnum=" << local_eventnum << endl;
		  // h2_leddown_z_vs_event->Fill(local_eventnum,z);    
		  h2_leddown_z_vs_event->Fill(eventnumber,z);   
		  h2_leddown_Eup_vs_Edown->Fill(Edown,Eup);   
		}
	}
				
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

    return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_BCAL_LED_time::erun(void) {
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


jerror_t JEventProcessor_BCAL_LED_time::fini(void) {
	// Called before program exit after event processing is finished.


return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
