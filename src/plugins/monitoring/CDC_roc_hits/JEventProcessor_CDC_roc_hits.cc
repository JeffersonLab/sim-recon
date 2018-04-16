// $Id$
//
//    File: JEventProcessor_CDC_roc_hits.cc
// Created: 18 May 2015
// Creator: Naomi Jarvis
//
// Plot channel number vs slot for any hit in each roc
// Useful to see which channels are very noisy, or dead
// Same purpose as occupancy, but organized by roc & slot for easier diagnostics
// Prints list of suspicious channels
//
// eg hd_root /raid12/gluex/rawdata2/Run003797/hd_rawdata_003797_000.evio -PPLUGINS=CDC_roc_hits -PEVIO:F125_FDC_WE=80 -PEVIO:F125_PI_EMULATION_MODE=0 -PEVIO:F125_CDC_WE=80 -o run3797.root
// 


#include <stdint.h>
#include <vector>
#include <stdio.h>

#include "JEventProcessor_CDC_roc_hits.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;




//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_roc_hits());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_roc_hits::JEventProcessor_CDC_roc_hits() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_roc_hits::~JEventProcessor_CDC_roc_hits() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_roc_hits::init(void) {

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!


  TSTART = 200;
  TSTOP = 1200;

  if (gPARMS) {
    gPARMS->SetDefaultParameter("CDC_ROC_HITS:TSTART",TSTART);
    gPARMS->SetDefaultParameter("CDC_ROC_HITS:TSTOP",TSTOP);
  }

  // create root folder for cdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("CDC_roc_hits")->cd();

  // book histograms
  cdc_nevents = new TH1I("cdc_nevents","CDC number of events",1, 0.5, 1.5);
  cdc_hits_roc25   = new TH2D("cdc_hits_roc25","CDC hits in ROC 25, channel vs slot;slot;channel",15,2.5,17.5,72,0.5,72.5);
  cdc_hits_roc26   = new TH2D("cdc_hits_roc26","CDC hits in ROC 26, channel vs slot;slot;channel",14,2.5,16.5,72,0.5,72.5);
  cdc_hits_roc27   = new TH2D("cdc_hits_roc27","CDC hits in ROC 27, channel vs slot;slot;channel",14,2.5,16.5,72,0.5,72.5);
  cdc_hits_roc28   = new TH2D("cdc_hits_roc28","CDC hits in ROC 28, channel vs slot;slot;channel",15,2.5,17.5,72,0.5,72.5);

  cdc_sumamp_roc25   = new TH1D("cdc_sumamp_roc25","CDC pulse peak amplitude summed for ROC 25;pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_roc26   = new TH1D("cdc_sumamp_roc26","CDC pulse peak amplitude summed for ROC 26;pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_roc27   = new TH1D("cdc_sumamp_roc27","CDC pulse peak amplitude summed for ROC 27;pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_roc28   = new TH1D("cdc_sumamp_roc28","CDC pulse peak amplitude summed for ROC 28;pulse peak (ADC units)",4096,0,4096);

  cdc_amp_roc25   = new TH2D("cdc_amp_roc25","CDC pulse peak amplitude in ROC 25;slot*100+channel;pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_roc26   = new TH2D("cdc_amp_roc26","CDC pulse peak amplitude in ROC 26;slot*100+channel;pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_roc27   = new TH2D("cdc_amp_roc27","CDC pulse peak amplitude in ROC 27;slot*100+channel;pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_roc28   = new TH2D("cdc_amp_roc28","CDC pulse peak amplitude in ROC 28;slot*100+channel;pulse peak (ADC units)",1600,200,1800,4096,0,4096);

  cdc_netamp_roc25   = new TH2D("cdc_netamp_roc25","CDC pulse peak amplitude minus pedestal in ROC 25;slot*100+channel;pulse peak minus pedestal (ADC units)",1600,200,1800,4096,0,4096);  
cdc_netamp_roc26   = new TH2D("cdc_netamp_roc26","CDC pulse peak amplitude minus pedestal in ROC 26;slot*100+channel;pulse peak minus pedestal (ADC units)",1600,200,1800,4096,0,4096);
cdc_netamp_roc27   = new TH2D("cdc_netamp_roc27","CDC pulse peak amplitude minus pedestal in ROC 27;slot*100+channel;pulse peak minus pedestal (ADC units)",1600,200,1800,4096,0,4096);
cdc_netamp_roc28   = new TH2D("cdc_netamp_roc28","CDC pulse peak amplitude minus pedestal in ROC 28;slot*100+channel;pulse peak minus pedestal (ADC units)",1600,200,1800,4096,0,4096);

  cdc_time_roc25   = new TH2D("cdc_time_roc25","CDC pulse time in ROC 25;slot*100+channel;time (0.8ns)",1600,200,1800,1000,0,2000);
  cdc_time_roc26   = new TH2D("cdc_time_roc26","CDC pulse time in ROC 26;slot*100+channel;time (0.8ns)",1600,200,1800,1000,0,2000);
  cdc_time_roc27   = new TH2D("cdc_time_roc27","CDC pulse time in ROC 27;slot*100+channel;time (0.8ns)",1600,200,1800,1000,0,2000);
  cdc_time_roc28   = new TH2D("cdc_time_roc28","CDC pulse time in ROC 28;slot*100+channel;time (0.8ns)",1600,200,1800,1000,0,2000);

  cdc_ped_roc25   = new TH2D("cdc_ped_roc25","CDC pedestal in ROC 25;slot*100+channel;pedestal",1600,200,1800,255,0,255);
  cdc_ped_roc26   = new TH2D("cdc_ped_roc26","CDC pedestal in ROC 26;slot*100+channel;pedestal",1600,200,1800,255,0,255);
  cdc_ped_roc27   = new TH2D("cdc_ped_roc27","CDC pedestal in ROC 27;slot*100+channel;pedestal",1600,200,1800,255,0,255);
  cdc_ped_roc28   = new TH2D("cdc_ped_roc28","CDC pedestal in ROC 28;slot*100+channel;pedestal",1600,200,1800,255,0,255);


  cdc_time_selected   = new TH1D("cdc_time_selected","CDC pulse time selected for amp_t histograms;time (0.8ns)",1000,0,2000);

 cdc_sumamp_t_roc25   = new TH1D("cdc_sumamp_t_roc25","CDC pulse peak amplitude summed for ROC 25 with time cut;pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_t_roc26   = new TH1D("cdc_sumamp_t_roc26","CDC pulse peak amplitude summed for ROC 26 with time cut; pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_t_roc27   = new TH1D("cdc_sumamp_t_roc27","CDC pulse peak amplitude summed for ROC 27 with time cut; pulse peak (ADC units)",4096,0,4096);
  cdc_sumamp_t_roc28   = new TH1D("cdc_sumamp_t_roc28","CDC pulse peak amplitude summed for ROC 28 with time cut; pulse peak (ADC units)",4096,0,4096);

  cdc_amp_t_roc25   = new TH2D("cdc_amp_t_roc25","CDC pulse peak amplitude in ROC 25 with time cut;slot*100+channel; pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_t_roc26   = new TH2D("cdc_amp_t_roc26","CDC pulse peak amplitude in ROC 26 with time cut;slot*100+channel; pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_t_roc27   = new TH2D("cdc_amp_t_roc27","CDC pulse peak amplitude in ROC 27 with time cut;slot*100+channel; pulse peak (ADC units)",1600,200,1800,4096,0,4096);
  cdc_amp_t_roc28   = new TH2D("cdc_amp_t_roc28","CDC pulse peak amplitude in ROC 28 with time cut;slot*100+channel; pulse peak (ADC units)",1600,200,1800,4096,0,4096);

  
  main->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!


  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

 
  uint32_t rocid;
  uint32_t slot;
  uint32_t channel;
  uint32_t amp;
  uint32_t ped;
  uint32_t time;
  uint32_t qf;
   
  int netamp;  // amp - pedestal

 // default scaling factors will be overridden by Df125Config if present
  uint16_t ASCALE = 8;   //amplitude
  uint16_t PSCALE = 1;   //ped


  // get raw data for cdc
  vector<const DCDCDigiHit*> digihits;
  eventLoop->Get(digihits);

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  if(digihits.size() > 0) cdc_nevents->Fill(1);


  for(uint32_t i=0; i<digihits.size(); i++) {

    const DCDCDigiHit *digihit = digihits[i]; 

    rocid = 0;
    slot = 0;
    channel = 0;
    amp = 0;
    ped = 0;
    time = 0;
    qf = 0;
    netamp = 0;

    const Df125PulseIntegral *pi = NULL; 
    digihit->GetSingle(pi);

    if (pi) {
      rocid = pi->rocid;
      slot = pi->slot;
      channel = pi->channel;
    }

    const Df125PulsePedestal *pp = NULL; 
    digihit->GetSingle(pp);

    if (pp) amp = pp->pulse_peak;

    const Df125CDCPulse *cp = NULL;
    digihit->GetSingle(cp);

    if (cp) {
      rocid = cp->rocid;
      slot = cp->slot;
      channel = cp->channel;
      amp = cp->first_max_amp;
      ped = cp->pedestal;
      time = cp->le_time;
      qf = cp->time_quality_bit;
      //      printf("cp amp %i\n",amp);

      const Df125Config *cf = NULL;

      cp->GetSingle(cf);
      if (cf) {
        ASCALE = 1<<cf->ABIT;
        PSCALE = 1<<cf->PBIT;
      }


    }


    netamp = ASCALE*amp - PSCALE*ped;
 


    if (!rocid) continue;

    if (rocid == 25) cdc_hits_roc25->Fill(slot,channel);
    if (rocid == 26) cdc_hits_roc26->Fill(slot,channel);
    if (rocid == 27) cdc_hits_roc27->Fill(slot,channel);
    if (rocid == 28) cdc_hits_roc28->Fill(slot,channel);

    if (rocid == 25) cdc_sumamp_roc25->Fill(amp);
    if (rocid == 26) cdc_sumamp_roc26->Fill(amp);
    if (rocid == 27) cdc_sumamp_roc27->Fill(amp);
    if (rocid == 28) cdc_sumamp_roc28->Fill(amp);

    if (rocid == 25) cdc_amp_roc25->Fill(slot*100+channel,amp);
    if (rocid == 26) cdc_amp_roc26->Fill(slot*100+channel,amp);
    if (rocid == 27) cdc_amp_roc27->Fill(slot*100+channel,amp);
    if (rocid == 28) cdc_amp_roc28->Fill(slot*100+channel,amp);

    if (rocid == 25) cdc_netamp_roc25->Fill(slot*100+channel,netamp);
    if (rocid == 26) cdc_netamp_roc26->Fill(slot*100+channel,netamp);
    if (rocid == 27) cdc_netamp_roc27->Fill(slot*100+channel,netamp);
    if (rocid == 28) cdc_netamp_roc28->Fill(slot*100+channel,netamp);

    if (rocid == 25) cdc_ped_roc25->Fill(slot*100+channel,ped);
    if (rocid == 26) cdc_ped_roc26->Fill(slot*100+channel,ped);
    if (rocid == 27) cdc_ped_roc27->Fill(slot*100+channel,ped);
    if (rocid == 28) cdc_ped_roc28->Fill(slot*100+channel,ped);

    if (!qf) {
      if (rocid == 25) cdc_time_roc25->Fill(slot*100+channel,time);
      if (rocid == 26) cdc_time_roc26->Fill(slot*100+channel,time);
      if (rocid == 27) cdc_time_roc27->Fill(slot*100+channel,time);
      if (rocid == 28) cdc_time_roc28->Fill(slot*100+channel,time);
    }

    if ((int)time>=TSTART && (int)time<=TSTOP) {

      cdc_time_selected->Fill(time);

      if (rocid == 25) cdc_sumamp_t_roc25->Fill(amp);
      if (rocid == 26) cdc_sumamp_t_roc26->Fill(amp);
      if (rocid == 27) cdc_sumamp_t_roc27->Fill(amp);
      if (rocid == 28) cdc_sumamp_t_roc28->Fill(amp);

      if (rocid == 25) cdc_amp_t_roc25->Fill(slot*100+channel,amp);
      if (rocid == 26) cdc_amp_t_roc26->Fill(slot*100+channel,amp);
      if (rocid == 27) cdc_amp_t_roc27->Fill(slot*100+channel,amp);
      if (rocid == 28) cdc_amp_t_roc28->Fill(slot*100+channel,amp);

    }

  }


  japp->RootUnLock(); //RELEASE ROOT LOCK!!


  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::fini(void) {
  // Called before program exit after event processing is finished.


  const float HOT = 3;
  const float COLD = 0.1;


  int maxslot;     // highest slot number used for roc in question
  int totalhits;   // total hits for roc
  int usedchans;   // number of channels used in roc
  int unusedchans; // number of unused channels in roc
  int meanhits = 1;    // mean number of hits per channel-in-use
  int chanhits;    // number of hits for current channel
  int nfound;      // number found in current search

  bool emptychan[15][72];
  bool hotchan[15][72];
  bool coldchan[15][72];
  
  bool unused;  // true if channel is not connected to a straw
  // nb there was one straw disconnected


  FILE *outfile;
  outfile = fopen("CDC_roc_hits.txt","w");

  if (outfile == NULL) printf("1 Cannot open file CDC_roc_hits.txt\n");
  if (!outfile) printf("2 Cannot open file CDC_roc_hits.txt\n");


  int roc;

  for (roc=25; roc<29; roc++) {

    for (int islot=3; islot<18; islot++) {
      for (int ichan=0; ichan<72; ichan++) {
        emptychan[islot-3][ichan] = 0;
        hotchan[islot-3][ichan] = 0;
        coldchan[islot-3][ichan] = 0;
      }
    }


    if (roc==25) {

      maxslot = 17;
      unusedchans = 14;
      totalhits = cdc_hits_roc25->GetEntries();

    } else if (roc==26) {

      maxslot = 16;
      unusedchans = 16;
      totalhits = cdc_hits_roc26->GetEntries();

    } else if (roc==27) {

      maxslot = 16;
      unusedchans = 16;
      totalhits = cdc_hits_roc27->GetEntries();

    } else if (roc==28) {

      maxslot = 17;
      unusedchans = 8;
      totalhits = cdc_hits_roc28->GetEntries();

    }

    usedchans = 72*(maxslot-4) - unusedchans;
    if (totalhits) meanhits = totalhits/usedchans;
    printf("\nROC %i: mean number of hits per channel is %i\n",roc,meanhits);
    if (outfile) fprintf(outfile,"\nROC %i: mean hits per channel is %i\n",roc,meanhits);



    for (int islot=3; islot<maxslot+1; islot++) {

      if (islot==11) continue;   // not used
      if (islot==12) continue;   // not used
 

      for (int ichan=0; ichan<72; ichan++) {

        unused = 0;

        if (roc==25) {
          if (islot==6 && ichan==17) unused = 1;
          if (islot==7 && ichan==11) unused = 1;
          if (islot==8 && ichan==35) unused = 1;
          if (islot==8 && ichan==37) unused = 1;
          if (islot==9 && ichan==37) unused = 1;
          if (islot==9 && ichan==59) unused = 1;
          if (islot==10 && ichan==61) unused = 1;
          if (islot==13 && ichan==13) unused = 1;
          if (islot==13 && ichan==37) unused = 1;
          if (islot==14 && ichan==35) unused = 1;
          if (islot==14 && ichan==59) unused = 1;
          if (islot==14 && ichan==25) unused = 1;  // straw W38 disconnected 17 Nov 2015
          if (islot==16 && ichan==11) unused = 1;
          if (islot==16 && ichan==37) unused = 1;
          if (islot==17 && ichan==59) unused = 1;

          chanhits = cdc_hits_roc25->GetBinContent(islot-2,ichan);

        } else if (roc==26) {

          if (islot==3 && ichan==59) unused = 1;
          if (islot==4 && ichan==7) unused = 1;
          if (islot==4 && ichan==11) unused = 1;
          if (islot==4 && ichan>47) unused = 1;
          if (islot==5 && ichan==13) unused = 1;
          if (islot==6 && ichan==11) unused = 1;
          if (islot==6 && ichan==61) unused = 1;
          if (islot==7 && ichan==13) unused = 1;
          if (islot==9 && ichan==33) unused = 1;
          if (islot==9 && ichan==35) unused = 1;
          if (islot==13 && ichan==37) unused = 1;
          if (islot==13 && ichan==61) unused = 1;
          if (islot==14 && ichan==61) unused = 1;
          if (islot==15 && ichan==61) unused = 1;
          if (islot==16 && ichan==13) unused = 1;
          if (islot==16 && ichan==35) unused = 1;

          chanhits = cdc_hits_roc26->GetBinContent(islot-2,ichan);

        } else if (roc==27) {

          if (islot==3 && ichan==37) unused = 1;
          if (islot==3 && ichan==61) unused = 1;
          if (islot==4 && ichan==11) unused = 1;
          if (islot==4 && ichan==35) unused = 1;
          if (islot==5 && ichan==13) unused = 1;
          if (islot==5 && ichan==35) unused = 1;
          if (islot==7 && ichan==17) unused = 1;
          if (islot==7 && ichan==41) unused = 1;
          if (islot==8 && ichan==13) unused = 1;  // straw K39 disconnected 6 Aug 2015
          if (islot==8 && ichan==61) unused = 1;
          if (islot==9 && ichan==11) unused = 1;
          if (islot==9 && ichan==35) unused = 1;
          if (islot==9 && ichan==37) unused = 1;
          if (islot==10 && ichan==37) unused = 1;
          if (islot==13 && ichan==35) unused = 1;
          if (islot==15 && ichan==59) unused = 1;
          if (islot==16 && ichan==11) unused = 1;

          chanhits = cdc_hits_roc27->GetBinContent(islot-2,ichan);

        } else if (roc==28) {

          if (islot==5 && ichan==37) unused = 1;
          if (islot==6 && ichan==65) unused = 1;
          if (islot==7 && ichan==13) unused = 1;
          if (islot==9 && ichan==35) unused = 1;
          if (islot==10 && ichan==37) unused = 1;
          if (islot==14 && ichan==31) unused = 1;
          if (islot==14 && ichan==37) unused = 1;
          if (islot==16 && ichan==59) unused = 1;
          if (islot==17 && ichan==61) unused = 1;

          chanhits = cdc_hits_roc28->GetBinContent(islot-2,ichan);

        }

        if (unused) continue;

        if (chanhits == 0) emptychan[islot-3][ichan] = 1;
        if (chanhits > HOT*meanhits) hotchan[islot-3][ichan] = 1;
        if (chanhits < COLD*meanhits) coldchan[islot-3][ichan] = 1;

      } //ichan
    }  //islot




    printf("Checking ROC %i for silent channels...\n",roc);
    if (outfile) fprintf(outfile,"Checking ROC %i for silent channels...\n",roc);


    nfound = 0;
    for (int islot=3; islot<maxslot; islot++) {
      for (int ichan=0; ichan<72; ichan++) {
        if (emptychan[islot-3][ichan]) {
          nfound++;
          printf("  roc %i slot %i channel %i\n",roc,islot,ichan);
          if (outfile) fprintf(outfile," roc %i slot %i channel %i\n",roc,islot,ichan);
        }
      }
    }

    if (!nfound) {
      printf("  none found \n");
      if (outfile) fprintf(outfile,"  none found \n");
    }

    

    printf("Checking ROC %i for cold channels (below %1.2f x mean hits/channel)...\n",roc,COLD);
    if (outfile) fprintf(outfile,"Checking ROC %i for cold channels (below %1.2f x mean hits/channel)...\n",roc,COLD);
    nfound = 0;

    for (int islot=3; islot<maxslot; islot++) {
      for (int ichan=0; ichan<72; ichan++) {
        if (coldchan[islot-3][ichan] && !emptychan[islot-3][ichan]) {
          nfound++;
          printf("  roc %i slot %i channel %i\n",roc,islot,ichan);
          if (outfile) fprintf(outfile," roc %i slot %i channel %i\n",roc,islot,ichan);
        }
      }
    }
    if (!nfound) {
      printf("  none found \n");
      if (outfile) fprintf(outfile,"  none found \n");
    }


    printf("Checking ROC %i for hot channels (over %1.2f x mean hits/channel)...\n",roc,HOT);
    if (outfile) fprintf(outfile,"Checking ROC %i for hot channels (over %1.2f x mean hits/channel)...\n",roc,HOT);
    nfound = 0;

    for (int islot=3; islot<maxslot; islot++) {
      for (int ichan=0; ichan<72; ichan++) {
        if (hotchan[islot-3][ichan]) {
          nfound++;
          printf("  roc %i slot %i channel %i\n",roc,islot,ichan);
          if (outfile) fprintf(outfile," roc %i slot %i channel %i\n",roc,islot,ichan);
        }
      }
    }
    if (!nfound) {
      printf("  none found \n");
      if (outfile) fprintf(outfile,"  none found \n");
    }



  } //roc


  if (outfile) fclose(outfile);


  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
