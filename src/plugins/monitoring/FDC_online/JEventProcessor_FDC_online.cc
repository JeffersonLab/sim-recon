// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>


#include "JEventProcessor_FDC_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "FDC/DFDCHit.h"
#include "FDC/DFDCWireDigiHit.h"
#include "FDC/DFDCCathodeDigiHit.h"
#include "DAQ/Df125PulsePedestal.h"
#include <DAQ/Df125PulseTime.h>

#include <TMath.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>

#define PI 3.14159265
#define TDC_v3BIN_SIZE       0.1118   // ns/LSB
//#define TDC_v3BIN_SIZE       1.   // ns/LSB

// root hist pointers
static TH1I *fdc_wire_occ[4][6];
static TH1I *fdc_cathode_occ[2][4][6];
static TH2I *fdc_cathode_time[2][4][6];
static TH2I *fdc_cathode_pulse_height[2][4][6];

static TH1I *fdc_num_events;

static TH2I *fdcow;
static TH2I *fdcos;

//------------------------------------------------------------------------------
// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_FDC_online());
  }
}


//------------------------------------------------------------------------------
JEventProcessor_FDC_online::JEventProcessor_FDC_online() {
}

//------------------------------------------------------------------------------
JEventProcessor_FDC_online::~JEventProcessor_FDC_online() {
}
//------------------------------------------------------------------------------
jerror_t JEventProcessor_FDC_online::init(void) {

  // soft threshold
  thresh=50; // threshold in f125 units (0-4096)

  wire_pitch=10.; // wire pitch in mm
  strip_pitch_u=5.007; // average upper strip pitch in mm
  strip_pitch_d=5.007; // average down strip pitch in mm
  strip_angle=15.*PI/180.; // strip angle now in rads
  cell_rot_step=60.*PI/180.; // cell rotation step now in rads

  // create root folder for fdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("FDC")->cd();

  // book hist
  fdc_num_events = new TH1I("fdc_num_events","FDC Number of events",1, 0.5, 1.5);
  fdcow   = new TH2I("fdcow","FDC wire occupancy by wire,gLayer",96,0.5,96.5,
		     24,0.5,24.5);
  fdcow->SetXTitle("wire number");
  fdcow->SetYTitle("plane");
  fdcos   = new TH2I("fdcos","FDC strip occupancy by strip,plane",192,0.5,
		     192.5,48,0.5,48.5);
  fdcos->SetXTitle("Strip number");
  fdcos->SetYTitle("plane");

  for (unsigned k=0;k<4;k++){
    char mypackagename[40];
    int package=k+1;
    sprintf(mypackagename,"Package_%d",package);
    gDirectory->mkdir(mypackagename)->cd();

    char hname[80],htitle[80];
    for (unsigned i=0;i<6;i++){
      int chamber=i+1;

      sprintf(hname,"fdc_pack%d_chamber%d_wire_occ",package,chamber);
      sprintf(htitle,"FDC wire occupancy for Package %d Chamber %d",package,
	      chamber);
      fdc_wire_occ[k][i]=new TH1I(hname,htitle,96,0.5,96.5);
      fdc_wire_occ[k][i]->SetXTitle("Wire number");

      sprintf(hname,"fdc_pack%d_chamber%d_upstream_cathode_occ",package,
	      chamber);
      sprintf(htitle,"FDC upstream cathode occupancy for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_occ[0][k][i]=new TH1I(hname,htitle,192,0.5,192.5);
      fdc_cathode_occ[0][k][i]->SetXTitle("Strip number");

       sprintf(hname,"fdc_pack%d_chamber%d_upstream_cathode_time",package,
	      chamber);
      sprintf(htitle,"FDC upstream cathode times for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_time[0][k][i]=new TH2I(hname,htitle,200,0,200,192,0.5,192.5);
      fdc_cathode_time[0][k][i]->SetXTitle("t");
      fdc_cathode_time[0][k][i]->SetYTitle("Strip number");
      
      sprintf(hname,"fdc_pack%d_chamber%d_upstream_cathode_pulse_height",package,
	      chamber);
      sprintf(htitle,"FDC upstream cathode pulse heights for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_pulse_height[0][k][i]=new TH2I(hname,htitle,100,0,1000,192,0.5,192.5);
      fdc_cathode_pulse_height[0][k][i]->SetXTitle("Pulse Height (ADC counts)");
      fdc_cathode_pulse_height[0][k][i]->SetYTitle("Strip number");


      sprintf(hname,"fdc_pack%d_chamber%d_downstream_cathode_occ",package,
	      chamber);
      sprintf(htitle,"FDC downstream cathode occupancy for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_occ[1][k][i]=new TH1I(hname,htitle,192,0.5,192.5);
      fdc_cathode_occ[1][k][i]->SetXTitle("Strip number");

      
      sprintf(hname,"fdc_pack%d_chamber%d_downstream_cathode_time",package,
	      chamber);
      sprintf(htitle,"FDC downstream cathode times for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_time[1][k][i]=new TH2I(hname,htitle,200,0,200,192,0.5,192.5);
      fdc_cathode_time[1][k][i]->SetXTitle("t");
      fdc_cathode_time[1][k][i]->SetYTitle("Strip number");

      sprintf(hname,"fdc_pack%d_chamber%d_downstream_cathode_pulse_height",package,
	      chamber);
      sprintf(htitle,"FDC downstream cathode pulse heights for Package %d Chamber %d",
	      package,chamber);
      fdc_cathode_pulse_height[1][k][i]=new TH2I(hname,htitle,100,0,1000,192,0.5,192.5);
      fdc_cathode_pulse_height[1][k][i]->SetXTitle("Pulse Height (ADC counts)");
      fdc_cathode_pulse_height[1][k][i]->SetYTitle("Strip number");
    }

    gDirectory->cd("../");
  } 
  
  // back to main dir
  main->cd();

  return NOERROR;
}


//------------------------------------------------------------------------------
jerror_t JEventProcessor_FDC_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}

//------------------------------------------------------------------------------
jerror_t JEventProcessor_FDC_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  for(int pack=0;pack<4;pack++){
    for(int cell=0;cell<6;cell++){
       for(int wire=0;wire<96;wire++){
          TDCnh[pack][cell][wire]=0;
          for(int hit=0;hit<20;hit++){
            TDCval[pack][cell][wire][hit]=-1000.;
          }
       }
       for(int ud=0;ud<2;ud++){
         for(int strip=0;strip<192;strip++){
            ADCnh[pack][cell][ud][strip]=0;
            for(int hit=0;hit<20;hit++){
              ADCmax[pack][cell][ud][strip][hit]=-1000.;
              ADCtime[pack][cell][ud][strip][hit]=-1000.;
            }
         } 
       } 
    }
  }

  // get anode digis
  vector<const DFDCWireDigiHit*>anodedigis;
  eventLoop->Get(anodedigis);

  // Get cathode digis
  vector<const DFDCCathodeDigiHit *>cathodedigis;
  eventLoop->Get(cathodedigis);

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
  
  if( (anodedigis.size()>0) || (cathodedigis.size()>0) )
	  fdc_num_events->Fill(1);

  for (unsigned int i=0;i<cathodedigis.size();i++){
    const DFDCCathodeDigiHit *digi=cathodedigis[i];

    int ud=0;
    if(digi->view==3)ud=1;
    int pack=digi->package-1;
    int cell=digi->chamber-1;
    int strip=digi->strip-1;

    const Df125PulsePedestal *pulseped;
    digi->GetSingle(pulseped);

    if (pulseped!=NULL){
      float peak=pulseped->pulse_peak;
      float ped=pulseped->pedestal;
      uint32_t p_pulse_num=pulseped->pulse_number;
      if(peak>ped+thresh){
	// cout<<" pack,cell,ud,strip,pn,peak="<<pack<<" "<<cell<<" "<<ud<<" "<<strip<<" "<<p_pulse_num<<" "<<peak-ped<<endl;
	ADCmax[pack][cell][ud][strip][ADCnh[pack][cell][ud][strip]]=peak-ped;

	const Df125PulseTime *pulsetime;
	digi->GetSingle(pulsetime);
	
	if(pulsetime!=NULL){ 
	  uint32_t t_pulse_num=pulsetime->pulse_number;
	  if(t_pulse_num==p_pulse_num){
	    //float time=pulsetime->time/64.;	    
	    float time=pulsetime->time/8;
	    ADCtime[pack][cell][ud][strip][ADCnh[pack][cell][ud][strip]]=time;
	  }
	}
	ADCnh[pack][cell][ud][strip]++;
      }
    }
  }
  
  for (unsigned int i=0;i<anodedigis.size();i++){
    const DFDCWireDigiHit *wdigi=anodedigis[i];
    int pack=wdigi->package-1;
    int cell=wdigi->chamber-1;
    int wire=wdigi->wire-1;
    int time=wdigi->time*TDC_v3BIN_SIZE/8.;
    TDCval[pack][cell][wire][TDCnh[pack][cell][wire]]=time;
    //     cout<<" pack,cell,wire="<<pack<<" "<<cell<<" "<<wire<<endl;
    TDCnh[pack][cell][wire]++;
  } 

  float val,val1,val2;
  for (int pack=0;pack<4;pack++){
     for (int cell=0;cell<6;cell++){
        for (int wire=0;wire<96;wire++){
           if(TDCval[pack][cell][wire][0]>0){
              fdc_wire_occ[pack][cell]->Fill(wire);

	      int plane=6*pack+cell+1;
	      fdcow->Fill(wire+1,plane);
           } // end if tdc hit
        } //end wire loop
     } //end cell loop
  } //end package loop

  for (int pack=0;pack<4;pack++){
     for (int cell=0;cell<6;cell++){
       for (int ud=0;ud<2;ud++){
	 for (int strip=0;strip<192;strip++){
	   val=ADCmax[pack][cell][ud][strip][0]; 
	   if(val>thresh){
	     val1=0; 
	     if(strip>2)val1=ADCmax[pack][cell][ud][strip-1][0]; 
	     val2=0; 
	     if(strip<191)val2=ADCmax[pack][cell][ud][strip+1][0]; 
	     if(val1>thresh||val2>thresh){
	       fdc_cathode_occ[ud][pack][cell]->Fill(strip);
	     
	       for (int nh=0;nh<ADCnh[pack][cell][ud][strip];nh++){
		 fdc_cathode_time[ud][pack][cell]->Fill(ADCtime[pack][cell][ud][strip][nh],strip);
		 fdc_cathode_pulse_height[ud][pack][cell]->Fill(ADCmax[pack][cell][ud][strip][nh],strip);

		 int plane=12*pack+2*cell+ud+1;
		 fdcos->Fill(strip+1,plane);
	       }
	     } //end if val>thresh
	   } //end if val1,2>thresh
	 } //end strip loop
       } //end ud loop
     } //end cell loop
  } //end package loop
  
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  return NOERROR;
}


//------------------------------------------------------------------------------
jerror_t JEventProcessor_FDC_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------------------------------------------------------------------
jerror_t JEventProcessor_FDC_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}

//------------------------------------------------------------------------------
