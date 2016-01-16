// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>
#include <limits>

#include "JEventProcessor_TOF_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "TOF/DTOFHit.h"
#include "TOF/DTOFDigiHit.h"
#include "TOF/DTOFTDCDigiHit.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPolyLine.h>


// root hist pointers
static TH1I *tofe;
static TH1I *toft;
static TH2I *tofo1;
static TH2I *tofo2;

static TH1I *tdcOccS;
static TH1I *tdcOccN;
static TH1I *tdcOccU;
static TH1I *tdcOccD;

static TH1I *adcOccS;
static TH1I *adcOccN;
static TH1I *adcOccU;
static TH1I *adcOccD;

static TH2I *planeHor;
static TH2I *planeVer;

static TH1I *tof_num_events;

static TH1I *histPed;
static TH1I *hTimeAdc;
static TH1I *hTimeTdc;

static TH2F *TOFPedestalsPlane0;
static TH2F *TOFPedestalsPlane1;
static TH2F *TOFSignalsRawPlane0;
static TH2F *TOFSignalsRawPlane1;
static TH2F *TOFTimesPlane0;
static TH2F *TOFTimesPlane1;

//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TOF_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_TOF_online::JEventProcessor_TOF_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TOF_online::~JEventProcessor_TOF_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TOF_online::init(void) {

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  // create root folder for tof and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("tof")->cd();


  // book hist
  tof_num_events = new TH1I("tof_num_events","TOF Number of events",1, 0.5, 1.5);

  tofe    = new TH1I("tofe","TOF energy in keV",100,0,5000);
  toft    = new TH1I("toft","TOF time in usec",100,0,500);
  tofo1   = new TH2I("tofo1","TOF occupancy plane 1 by bar,top/bottom",50,0,50,2,0,2);
  tofo2   = new TH2I("tofo2","TOF occupancy plane 2 by left/right,bar",2,0,2,50,0,50);

  tdcOccS = new TH1I("tdcOccS","TOF, TDC Occupancy",86,1,44);
  tdcOccN = new TH1I("tdcOccN","TOF, TDC Occupancy",86,1,44);
  tdcOccU = new TH1I("tdcOccU","TOF, TDC Occupancy",86,1,44);
  tdcOccD = new TH1I("tdcOccD","TOF, TDC Occupancy",86,1,44);

  adcOccS = new TH1I("adcOccS","TOF, fADC Occupancy",86,1,44);
  adcOccN = new TH1I("adcOccN","TOF, fADC Occupancy",86,1,44);
  adcOccU = new TH1I("adcOccU","TOF, fADC Occupancy",86,1,44);
  adcOccD = new TH1I("adcOccD","TOF, fADC Occupancy",86,1,44);

  histPed = new TH1I("histPed","TOF, Pedestals",50,190,210);

  hTimeAdc = new TH1I("hTimeAdc","TOF, fADC time",100,0,750);
  hTimeTdc = new TH1I("hTimeTdc","TOF, TDC time",100,0,750);

  planeHor = new TH2I("planeHor","TOF Upstream, Hit position, Horizontal Plane",84,-126,126,84,-126,126);
  planeVer = new TH2I("planeVer","TOF Upstream, Hit position, Vertical Plane",84,-126,126,84,-126,126);


  TOFPedestalsPlane0 = new TH2F("TOFPedestalsPlane0","TOF Pedestals Plane 0 all PMTs",
				100,50.,150., 88, 0., 88.);
  TOFPedestalsPlane1 = new TH2F("TOFPedestalsPlane1","TOF Pedestals Plane 1 all PMTs",
				100,50.,150., 88, 0., 88.);
  TOFSignalsRawPlane0 = new TH2F("TOFSignalsRawPlane0","TOF ADC Integral Plane 0 all PMTs",
				 300,0.,10000., 88, 0., 88.);
  TOFSignalsRawPlane1 = new TH2F("TOFSignalsRawPlane1","TOF ADC Integral Plane 1 all PMTs",
				 300,0.,10000., 88, 0., 88.);

  TOFTimesPlane0 = new TH2F("TOFTimesPlane0","TOF TDC times Plane 0 all PMTs",
			    800,0.,4000., 88, 0., 88.);
  TOFTimesPlane1 = new TH2F("TOFTimesPlane1","TOF TDC times Plane 1 all PMTs",
			    800,0.,4000., 88, 0., 88.);

  // back to main dir
  main->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TOF_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TOF_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  uint32_t E,t,pedestal;
  int plane,bar,end;
  int count_tdc = 0;
  int count_adc = 0;
  double tdc_time,pulse_time;
  
  double hit_north[45];
  double hit_south[45];
  double hit_up[45];
  double hit_down[45];
  double position, time, width;
  float integral;

  Float_t distY_Horz = -126; // Horizontal plane start counting from the Bottom to Top
  Float_t distX_Vert = -126; // Vertical plane start counting from the North to South

  memset(hit_north,0,sizeof(hit_north));
  memset(hit_south,0,sizeof(hit_south));
  memset(hit_up,0,sizeof(hit_up));
  memset(hit_down,0,sizeof(hit_down));

  // get data for tof
  vector<const DTOFHit*> dtofhits;
  eventLoop->Get(dtofhits);

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  if(dtofhits.size() > 0)
	  tof_num_events->Fill(1);

  for(unsigned int i=0; i<dtofhits.size(); i++) {
    const DTOFHit *dtofhit = dtofhits[i];

    plane = dtofhit->plane;
    bar   = dtofhit->bar;
    end   = dtofhit->end;
    E     = dtofhit->dE*1000000.;  // in MeV
    t     = dtofhit->t;         // in nanoseconds
    time     = dtofhit->t;         // in nanoseconds

    const DTOFDigiHit *digi = NULL;
    const DTOFTDCDigiHit *tdig = NULL;

    dtofhit->GetSingle(tdig);
    tdc_time = std::numeric_limits<double>::quiet_NaN(); //updated below if good
    pulse_time = std::numeric_limits<double>::quiet_NaN(); //updated below if good
    if(tdig){
      const uint32_t &Rtdc_time = tdig->time;
      if(Rtdc_time!=0xFFFF){ 
	tdc_time = 0.025*(double)Rtdc_time;
	if (plane){
	  TOFTimesPlane1->Fill((float)tdc_time, (float)bar-1+end*44);
	} else {
	  TOFTimesPlane0->Fill((float)tdc_time, (float)bar-1+end*44);
	}
      }
    }

    dtofhit->GetSingle(digi);
    pedestal = 0xFFFF;
    if (digi){
      const uint32_t &Rpedestal = digi->pedestal;
      const uint32_t &Rpulse_time = digi->pulse_time;
      if(Rpedestal!=0xFFFF && Rpulse_time!=0xFFFF){
        pedestal = Rpedestal;
	pulse_time = 0.0625*(double)Rpulse_time;
	integral = (float)digi->pulse_integral - (float)digi->pedestal*
	  (float)digi->nsamples_integral/(float)digi->nsamples_pedestal;
	integral /= (float)digi->nsamples_integral;
	integral *=10.;
	if (plane){
	  TOFPedestalsPlane1->Fill((float)pedestal,(float)bar-1+end*44);
	  TOFSignalsRawPlane1->Fill(integral, (float)bar-1+end*44);
	} else {
	  TOFPedestalsPlane0->Fill((float)pedestal,(float)bar-1+end*44);
	  TOFSignalsRawPlane0->Fill(integral, (float)bar-1+end*44);
	}
      }
    }
    

    if(dtofhits[i]->dE>0.0) {
      
      
      //      fill hist
      //      app->rootLock();
      
      tofe->Fill(E);
      toft->Fill(t);
      if(plane==0) {  // vertical, North is 0
        tofo1->Fill(bar,end);
      } else {
        tofo2->Fill(end,bar);
      }
                  
      //    app->rootUnlock();
      
      
    } // close if E>0

    if (dtofhit->has_fADC) hTimeAdc->Fill(pulse_time);
    if (dtofhit->has_TDC) hTimeTdc->Fill(tdc_time);

    switch(plane)
      {
      case 0:
	if(end == 1){
	  //cout << "Down : " << bar << endl;
	  if (dtofhit->has_fADC){
	    if (pedestal!=0xFFFF) histPed->Fill(pedestal);
	    adcOccD->Fill(bar);
	    count_adc++;
	  }
	  if (dtofhit->has_TDC){
	    tdcOccD->Fill(bar);
	    count_tdc++;
	    if ( ( (hit_down[bar]<=0) || (t < hit_down[bar]) ) && (bar!=22 || bar!=23) ){
	      hit_down[bar] = time;
	    }
	  }
	}
	else if(end == 0){
	  //cout << "Up : " << bar << endl;
	  if (dtofhit->has_fADC){
	    if (pedestal!=0xFFFF) histPed->Fill(pedestal);
	    adcOccU->Fill(bar);
	    count_adc++;
	  }
	  if (dtofhit->has_TDC){
	    tdcOccU->Fill(bar);
	    count_tdc++;
	    if ( ( (hit_up[bar]<=0) || (t < hit_up[bar]) ) && (bar!=22 || bar!=23) ){
	      hit_up[bar] = time;
	    }
	  }
	}
	break;
      case 1:
	if(end == 0){
	  //cout << "North : " << bar << endl;
	  if (dtofhit->has_fADC){
	    if (pedestal!=0xFFFF) histPed->Fill(pedestal);
	    adcOccN->Fill(bar);
	    count_adc++;
	  }
	  if (dtofhit->has_TDC){
	    tdcOccN->Fill(bar);
	    count_tdc++;
	    if ( ( (hit_north[bar]<=0) || (t < hit_north[bar]) ) && (bar!=22 || bar!=23) ){
	      hit_north[bar] = time;
	    }
	  }
	}
	else if(end == 1){
	  //cout << "South : " << bar << endl;
	  if (dtofhit->has_fADC){
	    if (pedestal!=0xFFFF) histPed->Fill(pedestal);
	    adcOccS->Fill(bar);
	    count_adc++;
	  }
	  if (dtofhit->has_TDC){
	    tdcOccS->Fill(bar);
	    count_tdc++;
	    if ( ( (hit_south[bar]<=0) || (t < hit_south[bar]) ) && (bar!=22 || bar!=23) ){
	      hit_south[bar] = time;
	    }
	  }
	}
	break;
      }
    
  } // close DTOFHit size


  for (int i=1; i<45; i++)
    {

      if ( i == 20 || i == 21 || i == 24 || i == 25 ){
	distY_Horz = distY_Horz + 1.5;
	distX_Vert = distX_Vert + 1.5;
	width = 3.0;
      }
      else{
	distY_Horz = distY_Horz + 3;
	distX_Vert = distX_Vert + 3;
	width = 6.0;
      }

      if( hit_south[i]>0 && hit_north[i]>0 )
	{
	  position = (15.2*(Float_t(hit_south[i] - hit_north[i])/2) );
	  distY_Horz = distY_Horz + (drand48()-0.5)*width;
	  if (position )planeHor->Fill(position,distY_Horz);
	}

      if( hit_up[i]>0 && hit_down[i]>0 )
	{
	  position = (15.2*(Float_t(hit_up[i] - hit_down[i])/2) );
	  distX_Vert = distX_Vert + (drand48()-0.5)*width;
	  planeVer->Fill(distX_Vert,position);
	}

      if ( i == 20 || i == 21 || i == 24 || i == 25 ){
	distY_Horz = distY_Horz + 1.5;
	distX_Vert = distX_Vert + 1.5;
      }
      else{
	distY_Horz = distY_Horz + 3;
	distX_Vert = distX_Vert + 3;
      }

    }

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TOF_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TOF_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
