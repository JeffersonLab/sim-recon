// $Id$
//
//    File: JEventProcessor_ST_online_multi.cc
// Created: Wed Feb 24 08:17:21 EST 2016
// Creator: mkamel (on Linux mkamel-NE56R 3.13.0-39-generic x86_64)
//

#include "JEventProcessor_ST_online_multi.h"
#include "TRIGGER/DTrigger.h"

using namespace jana;
//***************** Declare Two Dimensional Histograms*************
static TH2I *h2_st_adc_tdc_multi;
static TH2I *h2_st_adc_hit_multi;
static TH2I *h2_sector_adc_multip;
static TH2I *h2_sector_tdc_multip;
static TH2I *h2_sector_hit_multip;
static TH2I *h2_adc_unmatched;
static TH2I *h2_tdc_unmatched;
//***************** Declare One Dimensional Histograms****************
static TH1I *h1_adc_sec;
static TH1I *h1_tdc_sec;
static TH1I *h1_hit_sec;
static TH1I *st_num_events;
static TH1I *h1_tdc_noMatchTO_adc;
static TH1I *h1_adc_noMatchTO_tdc;
//***************** Declare Dynamic Arrays of  Histograms *************
TH1I** h1_ADC_multiplicity = new TH1I*[NCHANNELS];
TH1I** h1_TDC_multiplicity = new TH1I*[NCHANNELS];
TH1I** h1_hit_multiplicity = new TH1I*[NCHANNELS];
TH2I** h2_ADC_TDC_multiplicity = new TH2I*[NCHANNELS];

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_ST_online_multi());
}
} // "C"


//------------------
// JEventProcessor_ST_online_multi (Constructor)
//------------------
JEventProcessor_ST_online_multi::JEventProcessor_ST_online_multi()
{

}

//------------------
// ~JEventProcessor_ST_online_multi (Destructor)
//------------------
JEventProcessor_ST_online_multi::~JEventProcessor_ST_online_multi()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ST_online_multi::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

  //Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("st_multiplicity")->cd();
  //Define 1D histograms
  st_num_events = new TH1I("st_num_events","ST Number of events",1, 0.5, 1.5);
  h1_adc_sec = new TH1I("h1_adc_sec", "ST fADC250 DigiHit Occupancy; Channel Number; fADC250 Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);
  h1_tdc_sec = new TH1I("h1_tdc_sec", "ST TDC DigiHit Occupancy; Channel Number; TDC Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);
  h1_hit_sec = new TH1I("h1_hit_sec", "ST Hit Occupancy; Channel Number; Hit Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);

  h1_tdc_noMatchTO_adc = new TH1I("h1_tdc_noMatchTO_adc", "TDC not matched to ADC; Channel Number; TDC Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);
  h1_adc_noMatchTO_tdc = new TH1I("h1_adc_noMatchTO_tdc", "ADC not matched to TDC; Channel Number; ADC Counts", NCHANNELS, 0.5, NCHANNELS + 0.5);


  // 2D Multiplicity Histos
  h2_st_adc_tdc_multi = new TH2I("h2_st_adc_tdc_multi", "ST Total Multiplicity: TDC vs ADC; f1TDC Multiplicity; fADC250 Multiplicity", TDC_MULTI_BINS, TDC_MULTI_MIN + 0.5, TDC_MULTI_MAX + 0.5, ADC_MULTI_BINS, ADC_MULTI_MIN + 0.5, ADC_MULTI_MAX + 0.5);
  h2_st_adc_hit_multi = new TH2I("h2_st_adc_hit_multi", "ST Total Multiplicity: HIT vs ADC; Hit Multiplicity; fADC250 Multiplicity", TDC_MULTI_BINS, TDC_MULTI_MIN + 0.5, TDC_MULTI_MAX + 0.5, ADC_MULTI_BINS, ADC_MULTI_MIN + 0.5, ADC_MULTI_MAX + 0.5);
 
  h2_sector_adc_multip = new TH2I("h2_sector_adc_multip", "Sector Multiplicity: per event; Sector; ADC Multiplicity", NCHANNELS, 0.5, NCHANNELS + 0.5, 10, -0.5, 9.5);
  h2_sector_tdc_multip= new TH2I("h2_sector_tdc_multip", "Sector Multiplicity: per event; Sector; TDC Multiplicity", NCHANNELS, 0.5, NCHANNELS + 0.5, 10, -0.5, 9.5);
  h2_sector_hit_multip= new TH2I("h2_sector_hit_multip", "Sector Multiplicity: per event; Sector; Hit Multiplicity", NCHANNELS, 0.5, NCHANNELS + 0.5, 10, -0.5, 9.5);
  h2_adc_unmatched = new TH2I("h2_adc_unmatched", "Unmatched ADC to TDC: per event; Sector; ADC Unmatched Counts", NCHANNELS, 0.5, NCHANNELS + 0.5, 10, -0.5, 9.5);
  h2_tdc_unmatched = new TH2I("h2_tdc_unmatched", "Unmatched TDC to ADC: per event; Sector; TDC Unmatched Counts", NCHANNELS, 0.5, NCHANNELS + 0.5, 10, -0.5, 9.5);
  // 1d multiplicity per paddle per event
  for(int i = 0; i < NCHANNELS; i++)
    {
      h1_ADC_multiplicity[i] = new TH1I(Form("h1_multiplicity_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; fADC250 Multiplicity; fADC250 Counts ", i+1, 0+12*i, 12+12*i), 10, -0.5, 9.5);
      h1_TDC_multiplicity[i] =  new TH1I(Form("h1_TDC_multiplicity_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; TDC Multiplicity; TDC Counts ", i+1, 0+12*i, 12+12*i), 10, -0.5, 9.5);
      h1_hit_multiplicity[i] =  new TH1I(Form("h1_hit_multiplicity_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}]; Hit Multiplicity; Hit Counts ", i+1, 0+12*i, 12+12*i), 10, -0.5, 9.5); 
      // 2d multiplicity per paddle per event
        h2_ADC_TDC_multiplicity[i] = new TH2I (Form("h2_ADC_TDC_multiplicity_%i", i+1), Form("Channel %i, #phi #in [%i^{#circ}, %i^{#circ}];f1TDC Multiplicity; fADC250 Multiplicity", i+1, 0+12*i, 12+12*i), 10, -0.5, 9.5,  10, -0.5, 9.5);

    }

 
  // cd back to main directory
  main->cd();
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ST_online_multi::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ST_online_multi::evnt(JEventLoop *loop, uint64_t eventnumber)
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
  // Get the data objects first so we minimize the time we hold the ROOT mutex lock
  vector<const DSCDigiHit*> ADC_Digi_Hits;      // ST fADC250 DigiHits
  vector<const DSCTDCDigiHit*> TDC_Digi_Hits;   // ST f1TDC DigiHits
  vector<const DSCHit*> Factory_Hits;           // ST hits
  
  const DTrigger* locTrigger = NULL; 
  loop->GetSingle(locTrigger); 
  if(locTrigger->Get_L1FrontPanelTriggerBits() != 0)
    return NOERROR;

  loop->Get(ADC_Digi_Hits);
  loop->Get(TDC_Digi_Hits);
  loop->Get(Factory_Hits);
  //Get the size of each object
  int ADC_hits       = ADC_Digi_Hits.size();
  int TDC_hits       = TDC_Digi_Hits.size();
  int FAC_hits       = Factory_Hits.size();
 

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

 //reset the counters to zero
  memset(counter_adc, 0, sizeof(counter_adc));
  memset(counter_tdc, 0, sizeof(counter_tdc));
  memset(counter_hit, 0, sizeof(counter_hit));
  memset(counter_adc_2, 0, sizeof(counter_adc_2));
  memset(counter_tdc_2, 0, sizeof(counter_tdc_2));
  memset(counter_adc_unmatched,0, sizeof(counter_adc_unmatched));
  memset(counter_tdc_unmatched,0, sizeof(counter_tdc_unmatched));
//int hit_index;
  //Fill the number of events histograms
  if( (ADC_hits > 0) || (TDC_hits > 0) || (FAC_hits > 0) )
    st_num_events->Fill(1);
  //Fill 2D multiplicity histos
  h2_st_adc_tdc_multi->Fill(TDC_hits, ADC_hits);
  h2_st_adc_hit_multi->Fill(FAC_hits, ADC_hits);
  //=============================  ADC Hits ====================
  for(int i = 0; i < ADC_hits; i++)
    {
      int adc_sector  = ADC_Digi_Hits[i]->sector;
      adc_index   = adc_sector - 1;
      //Fill adc occupancy histo
      h1_adc_sec->Fill(adc_sector);
      //multiplicity per paddle per event
      counter_adc[adc_index] += 1;
      for(int k = 0; k < TDC_hits; k++)
	{ 
	  int tdc_sector = TDC_Digi_Hits[k]->sector;
	  if (adc_sector != tdc_sector)
	    {
	      counter_adc_unmatched[adc_index] +=1;
	      h1_tdc_noMatchTO_adc->Fill(tdc_sector);
	      h1_adc_noMatchTO_tdc->Fill(adc_sector);
	    }
	}
    }
  for(int j=0; j < NCHANNELS; j++)
    {
      if (counter_adc[j] != 0)
	{
	  h1_ADC_multiplicity[j]->Fill(counter_adc[j]);
	  h2_sector_adc_multip->Fill(j+1,counter_adc[j]);
	
      if ((counter_adc_unmatched[j] != 0) && (counter_adc_unmatched[j] <= 3))
	h2_adc_unmatched->Fill(j+1,counter_adc_unmatched[j]);
	}
    }

  //=============================  TDC Hits ====================
  // cout<< "**********************event number = "<< eventnumber <<"*****************************" << endl;  
  for(int k = 0; k < TDC_hits; k++)
    { 
      int tdc_sector = TDC_Digi_Hits[k]->sector;
      tdc_index = tdc_sector -1 ;
      //Fill tdc occupancy histo
      h1_tdc_sec->Fill(tdc_sector);
      //multiplicity per paddle per event
      counter_tdc[tdc_index] += 1;
      for(int i = 0; i < ADC_hits; i++)
	{ 
	  int adc_sector = ADC_Digi_Hits[i]->sector;
	  if (tdc_sector != adc_sector)
	    {
	      counter_tdc_unmatched[tdc_index] +=1;
	    }
	}

    }
  for(int j=0; j < NCHANNELS; j++)
    {
      if (counter_tdc[j] != 0)
	{
	  // cout << "counter_tdc[j] = " << counter_tdc[j] << endl; 
	  h1_TDC_multiplicity[j]->Fill(counter_tdc[j]);
	  h2_sector_tdc_multip->Fill(j+1,counter_tdc[j]);
	
      if ((counter_tdc_unmatched[j] != 0) && (counter_tdc_unmatched[j] <= 3))
	h2_tdc_unmatched->Fill(j+1,counter_tdc_unmatched[j]);
	}
    }
  //=========================== ST Factory Hits ====================
  for(int l = 0; l < FAC_hits; l++)
    {        
      int hit_sector = Factory_Hits[l]->sector;
      hit_index = hit_sector -1 ;
      //Fill tdc occupancy histo
      h1_hit_sec->Fill(hit_sector);
      //multiplicity per paddle per event
      counter_hit[hit_index] += 1;
    }    
  for(int j=0; j < NCHANNELS; j++)
    {
      if (counter_hit[j] != 0)
	{  
	  h1_hit_multiplicity[j]->Fill(counter_hit[j]);
	  h2_sector_hit_multip->Fill(j+1,counter_hit[j]);
	}
    }
  //==================== 2D multiplicity =============================
  for(int j=0; j < NCHANNELS; j++)
    {	  
      if ((counter_adc[j] != 0) && (counter_tdc[j] != 0))
	{  
	  h2_ADC_TDC_multiplicity[j]->Fill(counter_tdc[j],counter_adc[j]);
	}
    }
  
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ST_online_multi::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ST_online_multi::fini(void)
{
	// Called before program exit after event processing is finished.

 
	return NOERROR;
}

