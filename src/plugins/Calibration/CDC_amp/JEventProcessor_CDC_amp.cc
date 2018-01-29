// $Id$
//
//    File: JEventProcessor_CDC_amp.cc
// Created: Tue Sep  6 10:13:02 EDT 2016
// Creator: njarvis (on Linux egbert 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_CDC_amp.h"



// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_CDC_amp());
}
} // "C"


//------------------
// JEventProcessor_CDC_amp (Constructor)
//------------------
JEventProcessor_CDC_amp::JEventProcessor_CDC_amp()
{

}

//------------------
// ~JEventProcessor_CDC_amp (Destructor)
//------------------
JEventProcessor_CDC_amp::~JEventProcessor_CDC_amp()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_CDC_amp::init(void)
{
	// This is called once at program startup. 

  TDirectory *main = gDirectory;
  gDirectory->mkdir("CDC_amp")->cd();



  asum = new TH1I("asum","CDC amplitude (all hits);amplitude - pedestal",4096,0,4096);
  an = new TH2I("an","CDC amplitude vs n (all hits); n; amplitude - pedestal",3522,0,3522,4096,0,4096);

  atsum = new TH1I("atsum","CDC amplitude (hits on tracks);amplitude - pedestal",4096,0,4096);
  atn = new TH2I("atn","CDC amplitude vs n (hits on tracks); n; amplitude - pedestal",3522,0,3522,4096,0,4096);


  attsum = new TH1I("attsum","CDC amplitude (tracks, theta 85-95 deg, z 50-80cm);amplitude - pedestal",4096,0,4096);
  attn = new TH2I("attn","CDC amplitude vs n (tracks, theta 85-95 deg, z 50-80cm); n; amplitude - pedestal",3522,0,3522,4096,0,4096);

  atheta = new TH2D("atheta","CDC amplitude vs theta (hits on tracks); theta; amplitude - pedestal",180,0,180,4096,0,4096);


  qsum = new TH1D("qsum","charge (all hits);charge",1000,0,4e4);
  qn = new TH2D("qn","charge vs n (all hits); n; charge",3522,0,3522,1000,0,4e4);

  qtsum = new TH1D("qtsum","charge (hits on tracks);charge",1000,0,4e4);
  qtn = new TH2D("qtn","charge vs n (hits on tracks); n; charge",3522,0,3522,400,0,4e4);

  qttsum = new TH1D("qttsum","charge (tracks, theta 85-95 deg, z 50-80 cm);charge",1000,0,4e4);
  qttn = new TH2D("qttn","charge vs n (tracks, theta 85-95 deg, z 50-80 cm); n; charge",3522,0,3522,400,0,4e4);

  qtheta = new TH2D("qtheta","charge vs theta (hits on tracks); theta; charge",180,0,180,1000,0,4e4);



  main->cd();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_CDC_amp::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes

        if (runnumber<40000) ASCALE = 8;    // default for ASCALE before run 40,000 to be used if Df125config is not present

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_CDC_amp::evnt(JEventLoop *loop, uint64_t eventnumber)
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
	//  ... fill histograms or trees ...
	// japp->RootFillUnLock(this);

  int ring, straw, n;   // ring number, straw number within ring, straw number overall (1 to 3522)

  uint32_t amp,ped;     // dcdcdigihits raw quantities: time, pedestal, amplitude, quality factor, overflow count

  //scaling factors will be overridden by Df125Config if presqnt
  //  uint16_t ASCALE = 8;   //amplitude
  //  uint16_t PSCALE = 1;   //ped

  //add extra 0 at front to use offset[1] for ring 1  //used to calculate straw number n 
  int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};


  // select events with physics events, i.e., not LED and other front panel triggers
  const DTrigger* locTrigger = NULL; 
  loop->GetSingle(locTrigger); 
  if(locTrigger->Get_L1FrontPanelTriggerBits() != 0) 
    return NOERROR;

  // test whether this is simulated or real data (skip digihits for sim data)
  int SIMULATION;
  vector<const DMCThrown*> MCThrowns;
  loop->Get(MCThrowns);
  if (MCThrowns.empty()) SIMULATION = 0;
  if (!MCThrowns.empty()) SIMULATION = 1;


  vector<const DCDCHit*> hits;
  loop->Get(hits);

  const DCDCHit *hit = NULL;
  const DCDCDigiHit *digihit = NULL;
  const Df125CDCPulse *cp = NULL;
  const Df125Config *config = NULL;

  int netamp = 0;
  float scaledped;
  double charge;

  int used[3522] = {0};

  for (uint32_t i=0; i<hits.size(); i++) {

    hit = hits[i];
    netamp = 0;

    if (!SIMULATION) {

      digihit = NULL;
      hit->GetSingle(digihit);
      if (!digihit) continue;

      cp = NULL; 
      digihit->GetSingle(cp);
      if (!cp) continue; //no CDCPulseData (happens occasionally)

      if (cp->time_quality_bit) continue;
      if (cp->overflow_count) continue;

      config = NULL;
      cp->GetSingle(config);

      if (config) {  //defaults were set already in case config does not exist

        PSCALE = 1<<config->PBIT;
        ASCALE = 1<<config->ABIT;

      } 

      amp = cp->first_max_amp;
      ped = cp->pedestal;

      scaledped = ped*(float)PSCALE/(float)ASCALE;

      netamp = (int)amp - (int)scaledped;

    }

    charge = hit->q;

    ring     = (int)hit->ring;
    straw    = (int)hit->straw;

    n = straw_offset[ring] + straw;
 
    japp->RootFillLock(this); //ACQUIRE ROOT LOCK!!

    if (netamp>0) {
      asum->Fill(netamp);
      an->Fill(n,netamp);
    }

    if (charge>0) {
      qsum->Fill(charge);
      qn->Fill(n,charge);
    }

    japp->RootFillUnLock(this); //ACQUIRE ROOT LOCK!!


  }


  //--------tracks---------------------------
  

  vector<const DTrackTimeBased*> tracks;

  loop->Get(tracks);


  for (uint32_t i=0; i<tracks.size(); i++) {

    DVector3 mom = tracks[i]->momentum();
    double theta = mom.Theta();
    theta = 180.0*theta/3.14159;

    vector<DTrackFitter::pull_t> pulls = tracks[i]->pulls;

    for (uint32_t j=0; j<pulls.size(); j++) {

      if (pulls[j].cdc_hit == NULL) continue;

      hit = NULL;
      pulls[j].cdc_hit->GetSingle(hit);

      netamp = 0;  

      if (!SIMULATION) {

        digihit = NULL;
        hit->GetSingle(digihit);
        if (!digihit) continue;

        cp = NULL; 
        digihit->GetSingle(cp);
        if (!cp) continue; //no CDCPulseData (happens occasionally)

        if (cp->time_quality_bit) continue;
        if (cp->overflow_count) continue;

        amp = cp->first_max_amp;
        ped = cp->pedestal;

        scaledped = ped*(float)PSCALE/(float)ASCALE;

        netamp = (int)amp - (int)scaledped;

      }

      ring     = (int)hit->ring;
      straw    = (int)hit->straw;

      n = straw_offset[ring] + straw;
 
      charge = hit->q;

      japp->RootFillLock(this); //ACQUIRE ROOT LOCK!!

      if (!used[n]) {

        used[n] = 1;

        if (netamp > 0) {
          atsum->Fill(netamp);
          atn->Fill(n,netamp);
          atheta->Fill(theta,netamp);
        }

        if (charge > 0) {
          qtsum->Fill(charge);
          qtn->Fill(n,charge);
          qtheta->Fill(theta,charge);
        }

        double z = pulls[j].z;

        if ((z>50) && (z<80) && (theta>85) && (theta<95)) {

          if (netamp > 0) {
            attsum->Fill(netamp);
            attn->Fill(n,netamp);
          }

          if (charge > 0) {
            qttsum->Fill(charge);
            qttn->Fill(n,charge);
          }
        }

      }

      japp->RootFillUnLock(this); 

    }
  }



	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_CDC_amp::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_CDC_amp::fini(void)
{
	// Called before program exit after event processing is finished.

  if (!asum->GetEntries()) delete asum;
  if (!an->GetEntries()) delete an;

  if (!atsum->GetEntries()) delete atsum;
  if (!atn->GetEntries()) delete atn;

  if (!attsum->GetEntries()) delete attsum;
  if (!attn->GetEntries()) delete attn;

  if (!atheta->GetEntries()) delete atheta;


	return NOERROR;
}
