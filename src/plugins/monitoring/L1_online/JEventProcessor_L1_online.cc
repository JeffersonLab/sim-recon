// $Id$
//
//    File: JEventProcessor_L1_online.cc
// Created: Fri Mar 20 16:32:04 EDT 2015
//
#include <iostream>
#include <sstream>
#include "JEventProcessor_L1_online.h"

using namespace std;
using namespace jana;

#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <TRIGGER/DL1Trigger.h>

#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulsePedestal.h>

#include <FCAL/DFCALDigiHit.h>
#include <BCAL/DBCALDigiHit.h>

#include "START_COUNTER/DSCHit.h"

#include <TAGGER/DTAGHDigiHit.h>


// root hist pointers
static TH1I *htrig_bit;
static TH1I *htrig_bit_fp;

static TH1I *htrig_type;

static TH1F *hrate_gtp[8];
// static TH1F *hrate_fp[16];

static TH1I *hfcal_time1;
static TH1I *hfcal_time6;
static TH1I *hfcal_time7;

static TH1I *hbcal_time1;
static TH1I *hbcal_time3;
static TH1I *hbcal_time6;

static TH1I *htagh_time2;
static TH1I *htagh_occup2;

static TH1I *hfcal_en1;
static TH1I *hfcal_en6;
static TH1I *hfcal_en7;

static TH1I *hfcal_en2;

static TH1I *hbcal_en1;
static TH1I *hbcal_en3;
static TH1I *hbcal_en6;

static TH2I *hfcal_bcal_en1;
static TH2I *hfcal_bcal_en6;

static TH2I *hfcal_bcal_en1_3;

static TH2I *hfcal_bcal_type1;
static TH2I *hfcal_bcal_type3;
static TH2I *hfcal_bcal_type5;
static TH2I *hfcal_bcal_type7;

static TH1I *hst_hit2;
static TH1I *hst_hit7;

static TProfile *hrate_bit[8];


// static TH2F *hfcal_row_col;

//-------------------------
// Routine used to create our JEventProcessor
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_L1_online());

    }
} // "C"


//------------------
// JEventProcessor_L1_online (Constructor)
//------------------
JEventProcessor_L1_online::JEventProcessor_L1_online()
{

}

//------------------
// ~JEventProcessor_L1_online (Destructor)
//------------------
JEventProcessor_L1_online::~JEventProcessor_L1_online()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_L1_online::init(void)
{

    
  // FCAL constants - will be retrieved from the RCDB

    fcal_cell_thr  =  64;
    bcal_cell_thr  =  20;

  // Default trigger masks corresponding to 2 rings masked 
    fcal_row_mask_min = 26;
    fcal_row_mask_max = 32;
    fcal_col_mask_min = 26;
    fcal_col_mask_max = 32;

    const int bin_time = 100;
    const int bin_fcal = 100;
    const int bin_bcal = 200;
    


    // create root folder for pspair and cd to it, store main dir
    TDirectory *mainDir = gDirectory;
    TDirectory *L1Dir   = gDirectory->mkdir("L1");
    L1Dir->cd();

    // book hists
    htrig_bit = new TH1I("trig_bit","trig_bit", 32, 0.5, 32.5);
    htrig_bit_fp = new TH1I("trig_bit_fp","trig_bit_fp", 32, 0.5, 32.5);

    htrig_type = new TH1I("trig_type","trig_type", 257, -0.5, 256.5);

    for(Int_t ii = 0; ii < 8; ii++){
      
      //      char title[30];
      //      sprintf(title,"live_gtp_%d",ii);
      //      hlive_gtp[ii] = new TH1F(title,title,100,0.,100);
      
      char title[30];
      sprintf(title,"rate_gtp_%d",ii);
      hrate_gtp[ii] = new TH1F(title,title,500,0.,50.);      
      
      //      char title1[30];
      //      sprintf(title1,"rate_fp_%d",ii);
      //      hrate_fp[ii] = new TH1F(title1,title1,800,0.,400);     


      char title2[30];
      sprintf(title2,"rate_time_%d",ii);
      hrate_bit[ii] =  new TProfile(title2,title2, 501 ,0.,501,"s");

    }



    hfcal_time1  = new TH1I("fcal_time1","fcal_time1", bin_time, -0.5, 99.5);
    hfcal_time6  = new TH1I("fcal_time6","fcal_time6", bin_time, -0.5, 99.5); 
    hfcal_time7  = new TH1I("fcal_time7","fcal_time7", bin_time, -0.5, 99.5); 

    hbcal_time1  = new TH1I("bcal_time1","bcal_time1", bin_time, -0.5, 99.5);
    hbcal_time3  = new TH1I("bcal_time3","bcal_time3", bin_time, -0.5, 99.5); 
    hbcal_time6  = new TH1I("bcal_time6","bcal_time6", bin_time, -0.5, 99.5); 

    htagh_time2  = new TH1I("tagh_time2","tagh_time2", bin_time, -0.5, 99.5); 

    htagh_occup2 = new TH1I("tagh_occup2","tagh_occup2", 250, -0.5, 249.5); 

    hfcal_en1   =  new TH1I("fcal_en1", "fcal_en1", bin_fcal,  0., 20000.);
    hfcal_en2   =  new TH1I("fcal_en2", "fcal_en2", bin_fcal,  0., 20000.);

    hfcal_en6   =  new TH1I("fcal_en6", "fcal_en6", bin_fcal,  0., 20000.); 
    hfcal_en7   =  new TH1I("fcal_en7", "fcal_en7", bin_fcal,  0., 20000.); 

    hbcal_en1   =  new TH1I("bcal_en1", "bcal_en1", bin_bcal,  0., 50000.);
    hbcal_en3   =  new TH1I("bcal_en3", "bcal_en3", bin_bcal,  0., 50000.); 
    hbcal_en6   =  new TH1I("bcal_en6", "bcal_en6", bin_bcal,  0., 50000.); 

    hfcal_bcal_en1  = new TH2I("fcal_bcal_en1","fcal_bcal_en1", bin_fcal, 0., 20000, bin_bcal, 0., 50000);
    hfcal_bcal_en6  = new TH2I("fcal_bcal_en6","fcal_bcal_en6", bin_fcal, 0., 20000, bin_bcal, 0., 50000);

    hfcal_bcal_en1_3 = new TH2I("fcal_bcal_en_type1_3","fcal_bcal_en1_3", 100, 0., 20000, 200, 0., 50000);

    hfcal_bcal_type1 = new TH2I("fcal_bcal_type1","fcal_bcal_type1", bin_fcal, 0., 20000, bin_bcal, 0., 50000);
    hfcal_bcal_type3 = new TH2I("fcal_bcal_type3","fcal_bcal_type3", bin_fcal, 0., 20000, bin_bcal, 0., 50000);
    hfcal_bcal_type5 = new TH2I("fcal_bcal_type5","fcal_bcal_type5", bin_fcal, 0., 20000, bin_bcal, 0., 50000);
    hfcal_bcal_type7 = new TH2I("fcal_bcal_type7","fcal_bcal_type7", bin_fcal, 0., 20000, bin_bcal, 0., 50000);

    hst_hit7    =  new TH1I("st_hit7", "st_hit7", 100, -0.5, 99.5); 
    hst_hit2    =  new TH1I("st_hit2", "st_hit2", 100, -0.5, 99.5); 

    //    hfcal_row_col  = new TH2F("fcal_row_col","fcal_row_col", 60, 0.5, 60.5, 60, 0.5, 60.5);

    mainDir->cd();

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_L1_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  
    run_number = runnumber;
  
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_L1_online::evnt(JEventLoop *loop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop->Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.
  
     
  
      cout << " Run_number = " << run_number << endl;
      
      if( run_number < 11127 ){
	fcal_row_mask_min = 24;
	fcal_row_mask_max = 34;
	fcal_col_mask_min = 24;
	fcal_col_mask_max = 34;
      }

 


      vector<const DL1Trigger*>   l1trig; 

      vector<const DBCALDigiHit*> bcal_hits;
      vector<const DFCALDigiHit*> fcal_hits;      

      vector<const DSCHit*> st_hits;       

      vector<const DTAGHDigiHit*> tagh_hits;


      loop->Get(l1trig);
      loop->Get(bcal_hits);
      loop->Get(fcal_hits);

      loop->Get(st_hits);

      loop->Get(tagh_hits);

      memset(trig_bit,0,sizeof(trig_bit));
      memset(trig_bit_fp,0,sizeof(trig_bit_fp));


      // FILL HISTOGRAMS
      japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

      if( l1trig.size() > 0 ){
	for(unsigned int bit = 0; bit < 32; bit++){
	  trig_bit[bit + 1] = (l1trig[0]->trig_mask & (1 << bit)) ? 1 : 0;       
	  if(trig_bit[bit + 1] == 1){
	    htrig_bit->Fill(Float_t(bit+1)); 
	  }
	} 

	htrig_type->Fill(Float_t(l1trig[0]->trig_mask));
      
	for(unsigned int bit = 0; bit < 32; bit++){
	  trig_bit_fp[bit + 1] = (l1trig[0]->fp_trig_mask & (1 << bit)) ? 1 : 0;         
	  if(trig_bit_fp[bit + 1] == 1) htrig_bit_fp->Fill(Float_t(bit+1)); 
	} 
     
	// Sync Events
	if( l1trig[0]->gtp_rate.size() > 0 ){

	  for(unsigned int ii = 0; ii < 8; ii++){
	    hrate_gtp[ii]->Fill(Float_t(l1trig[0]->gtp_rate[ii]/1000.));
	    if( (l1trig[0]->nsync) > 0 && (l1trig[0]->nsync) < 500){
	      hrate_bit[ii]->Fill(Float_t(l1trig[0]->nsync),Float_t(l1trig[0]->gtp_rate[ii]/1000.));
	    }
	  }

	  //	  for(unsigned int ii = 0; ii < 16; ii++){
	  //	    hrate_fp[ii]->Fill(Float_t(l1trig[0]->fp_rate[ii]));
	  //	  }
	  

	}

      }      



      //-------------------   FCAL  ----------------------------
	
      int fcal_debug   = 0;      
      int fcal_tot_en  = 0;


      if(fcal_debug){
	cout << " Number of FCAL hits = " <<  fcal_hits.size() << endl;
      }
	

      for(unsigned int jj = 0; jj < fcal_hits.size(); jj++){
	
	const DFCALDigiHit *fcal_hit = fcal_hits[jj];	  
	
	const Df250PulsePedestal *pulsepedestal;
	const Df250PulseIntegral *pulseintegral;
	const Df250WindowRawData *windorawdata;
	
	fcal_hit->GetSingle(pulseintegral);	  
	pulseintegral->GetSingle(windorawdata);
	
	
	uint32_t adc_time   =    0;
	Int_t pulse_int     =   -10;
	Int_t pulse_peak    =   -10;
	

	adc_time = (fcal_hit->pulse_time & 0x7FC0) >> 6;

	pulse_int = fcal_hit->pulse_integral - fcal_hit->nsamples_integral*100;
	

	//Int_t fcal_raw_int  =  0;
	//Int_t fcal_peak = -10;
	//Int_t fcal_raw_time = -10;
	
	fcal_hit->GetSingle(pulsepedestal); 
	
	if(pulsepedestal){
	  
	  pulse_peak = pulsepedestal->pulse_peak - 100;
	  
	  if(fcal_debug){
	    cout << " Pulse peak      = " << pulse_peak << endl;
	    cout << " Pulse integral  = " << pulse_int << endl;
	    cout << " Pulse time      = " << adc_time  << endl;
	  }
	}
	
	int fcal_cell_used = 1;
	
	int row = fcal_hit->row; 
	int col = fcal_hit->column; 
	
	
	if( (row >= fcal_row_mask_min && row <= fcal_row_mask_max) 
	    && (col >= fcal_col_mask_min && col <= fcal_col_mask_max)) fcal_cell_used = 0;

	//	if(trig_bit[7] == 1)
	//	  hfcal_row_col->Fill(Float_t(fcal_hit->row),Float_t(fcal_hit->column));


	if(fcal_cell_used == 0) continue;
	

	if(pulse_peak > fcal_cell_thr){
	  fcal_tot_en += pulse_int;

	  if(trig_bit[1] == 1) hfcal_time1->Fill(Float_t(adc_time));
	  if(trig_bit[6] == 1) hfcal_time6->Fill(Float_t(adc_time));
	  if(trig_bit[7] == 1) hfcal_time7->Fill(Float_t(adc_time));
	  
	}	

      }



      //-------------------   BCAL  ----------------------------
	
      int bcal_debug  = 0;
      int bcal_tot_en = 0;
      
      if(bcal_debug){
	cout << " Number of BCAL hits = " <<  bcal_hits.size() << endl;
      }
      
	
      for(unsigned int jj = 0; jj < bcal_hits.size(); jj++){
	
	const DBCALDigiHit *bcal_hit = bcal_hits[jj];	  
	  
	const Df250PulsePedestal *pulsepedestal;
	const Df250PulseIntegral *pulseintegral;
	const Df250WindowRawData *windorawdata;
	
	bcal_hit->GetSingle(pulseintegral);	  
	pulseintegral->GetSingle(windorawdata);


	uint32_t adc_time   =    0;
	Int_t pulse_int     =   -10;
	Int_t pulse_peak    =   -10;

	adc_time = (bcal_hit->pulse_time & 0x7FC0) >> 6;

	pulse_int = bcal_hit->pulse_integral - bcal_hit->nsamples_integral*100;

	//Int_t bcal_raw_int   =   0;
	//Int_t bcal_peak      =  -10;
	//Int_t bcal_raw_time  =  -10;

	bcal_hit->GetSingle(pulsepedestal); 
	  
	if(pulsepedestal){
	  pulse_peak = pulsepedestal->pulse_peak - 100;
	  
	  if(bcal_debug){
	    cout << " Pulse peak      = " << pulse_peak << endl;
	    cout << " Pulse integral  = " << pulse_int << endl;
	    cout << " Pulse time      = " << adc_time  << endl;
	  }
	}
	
	if(pulse_peak > bcal_cell_thr){
	  if(trig_bit[1] == 1) hbcal_time1->Fill(Float_t(adc_time));
	  if(trig_bit[3] == 1) hbcal_time3->Fill(Float_t(adc_time));
	  if(trig_bit[6] == 1) hbcal_time6->Fill(Float_t(adc_time));

	  bcal_tot_en += pulse_int;
	}
      }


      if(trig_bit[1] == 1){
	hfcal_en1->Fill(Float_t(fcal_tot_en));
	hbcal_en1->Fill(Float_t(bcal_tot_en));
	hfcal_bcal_en1->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
      }

      if(trig_bit[6] == 1){
	hfcal_en6->Fill(Float_t(fcal_tot_en));
	hbcal_en6->Fill(Float_t(bcal_tot_en));
	hfcal_bcal_en6->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
      }

      if( (trig_bit[1] == 1) &&  (trig_bit[3] == 1) ){
      	hfcal_bcal_en1_3->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
      }

      if(trig_bit[7] == 1){
	hfcal_en7->Fill(Float_t(fcal_tot_en));
	hst_hit7->Fill(Float_t(st_hits.size()));
      }

      if(trig_bit[3] == 1){
	hbcal_en3->Fill(Float_t(bcal_tot_en));
      }   


      if(trig_bit[2] == 1){
	hfcal_en2->Fill(Float_t(fcal_tot_en));
      }   
      
      if( l1trig.size() > 0 ){
	if(l1trig[0]->trig_mask == 1)  hfcal_bcal_type1->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
	if(l1trig[0]->trig_mask == 3)  hfcal_bcal_type3->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
	if(l1trig[0]->trig_mask == 5)  hfcal_bcal_type5->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
	if(l1trig[0]->trig_mask == 7)  hfcal_bcal_type7->Fill(Float_t(fcal_tot_en), Float_t(bcal_tot_en));
      }


      // TAGH & ST

      for(unsigned int kk = 0; kk < tagh_hits.size(); kk++){
	const DTAGHDigiHit *tagh_hit = tagh_hits[kk];   
	const Df250PulseIntegral *pulseintegral;
        
	tagh_hit->GetSingle(pulseintegral);
	const Df250WindowRawData *windorawdata;
	pulseintegral->GetSingle(windorawdata);
        
	int counter_id = tagh_hit->counter_id;
	

	
	uint32_t adc_time = (tagh_hit->pulse_time & 0x7FC0) >> 6;
	
	if(trig_bit[2] == 1){
	  
	  if(counter_id < 100)
	    htagh_time2->Fill(Float_t(adc_time));

	  if( (adc_time > 5) && (adc_time < 15))
	    htagh_occup2->Fill(Float_t(counter_id));
	}
	
      }
      
      if(trig_bit[2] == 1){
	hst_hit2->Fill(Float_t(st_hits.size()));
      }



      japp->RootFillUnLock(this);   //RELEASE ROOT FILL LOCK
      
      return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_L1_online::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_L1_online::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

