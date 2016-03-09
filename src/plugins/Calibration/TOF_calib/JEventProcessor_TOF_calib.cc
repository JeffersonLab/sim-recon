// $Id$
//
//    File: JEventProcessor_TOF_calib.cc
// Created: Fri Mar 13 10:37:27 EDT 2015
// Creator: zihlmann (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_TOF_calib.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TOF_calib());
  }
} // "C"


//------------------
// JEventProcessor_TOF_calib (Constructor)
//------------------
JEventProcessor_TOF_calib::JEventProcessor_TOF_calib()
{
  
}

//------------------
// ~JEventProcessor_TOF_calib (Destructor)
//------------------
JEventProcessor_TOF_calib::~JEventProcessor_TOF_calib()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TOF_calib::init(void)
{
  // This is called once at program startup. If you are creating
  // and filling historgrams in this plugin, you should lock the
  // ROOT mutex like this:
  //
  // japp->RootWriteLock();
  //  ... fill historgrams or trees ...
  // japp->RootUnLock();
  //
  
  cout<<"INITIALIZE VARIABLES "<<flush<<endl;

  pthread_mutex_init(&mutex, NULL);

  //BINTDC_2_TIME = 0.025;
  BINTDC_2_TIME = 0.0234375;
  BINADC_2_TIME = 0.0625; // is 4ns/64

  TDCTLOC = 419.;
  ADCTLOC = 190.;

  ADCTimeCut = 50.;
  TDCTimeCut = 60.;

  first = 1;
  MakeHistograms();

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TOF_calib::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes

  RunNumber = runnumber;



  if (first){
    cout<<"SETUP HISTOGRAMS AND TREE FOR RUN IN brun()"<<RunNumber<<flush<<endl;
    MakeHistograms();
  }
  map<string,double> tdcshift;
  if (!eventLoop->GetCalib("/TOF/tdc_shift", tdcshift)){
    TOF_TDC_SHIFT = tdcshift["TOF_TDC_SHIFT"];
  }

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TOF_calib::evnt(JEventLoop *loop, uint64_t eventnumber)
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


  //cout<<"CALL EVENT ROUTINE!!!! "<<eventnumber<<endl;
  if (first){
    MakeHistograms();
  }
  

  // Get First Trigger Type
  const DL1Trigger *trig_words = NULL;
  uint32_t trig_mask, fp_trig_mask;
  try {
    loop->GetSingle(trig_words);
  } catch(...) {};
  if (trig_words) {
    trig_mask = trig_words->trig_mask;
    fp_trig_mask = trig_words->fp_trig_mask;
  }
  else {
    trig_mask = 0;
    fp_trig_mask = 0;
  }
  
  
  /* fp_trig_mask & 0x100 - upstream LED
     fp_trig_mask & 0x200 - downstream LED
     trig_mask & 0x1 - cosmic trigger*/
  
  if (fp_trig_mask){ // this is a front pannel trigger like LED
    return NOERROR;    
  }
  if (trig_mask>7){ // this is not a BCAL/FCAL trigger
    return NOERROR;    
  }
  
  vector<const DCODAEventInfo*> locCODAEventInfo;
  loop->Get(locCODAEventInfo);
  
  if (locCODAEventInfo.size() == 0){
    return NOERROR;
  }

  uint64_t TriggerTime = locCODAEventInfo[0]->avg_timestamp;
  int TriggerBIT = TriggerTime%6;
  float TimingShift = TOF_TDC_SHIFT - (float)TriggerBIT;
  if(TimingShift <= 0) { 
    TimingShift += 6.;
  } 
  
  TimingShift *= 4. ;

  vector<const DTOFDigiHit*> ADCHits, ADCHitsLeft[2],ADCHitsRight[2];
  vector<const DTOFTDCDigiHit*> TDCHits, TDCHitsLeft[2], TDCHitsRight[2];
  vector<int> ADCLeftOverFlow[2],ADCRightOverFlow[2];

  loop->Get(ADCHits);
  loop->Get(TDCHits);

  // loop over DTOFDigiHit: this are ADC hits 
  // sort them into ADCLeft and ADCRight hits
  // only keep hits within the time-peak
  // also keep the hodoscope planes separate 
  for (unsigned int k=0; k<ADCHits.size(); k++){
    const DTOFDigiHit *hit = ADCHits[k];
    int plane = hit->plane;
    int end = hit->end;
    float time = (float)hit->pulse_time * BINADC_2_TIME;
    int val = (hit->pulse_time & 0x3F);

    if (!val){
      continue;
    }

    TOFADCtime->Fill(time);
    if (fabsf(time-ADCTLOC)<ADCTimeCut){
      // test for overflow if raw data available
      vector <const Df250PulseIntegral*> PulseIntegral;
      hit->Get(PulseIntegral);
      vector <const Df250WindowRawData*> WRawData;
      PulseIntegral[0]->Get(WRawData);
      int overflow = 0;
      if ( WRawData.size() > 0) {
	for (int n=0;n<100;n++){
	  if (WRawData[0]->samples[n] & 0x1000){
	    overflow++;
	  }
	}
      } else {
	//overflow = PulseIntegral[0]->quality_factor;
      }
      
      if (end){
	ADCHitsRight[plane].push_back(hit);
	ADCRightOverFlow[plane].push_back(overflow);
      } else {
	ADCHitsLeft[plane].push_back(hit);
	ADCLeftOverFlow[plane].push_back(overflow);
      }
    }
  }

  if (ADCLeftOverFlow[0].size() != ADCHitsLeft[0].size()){
    cout<<"Error vector size missmatch!"<<endl;
  }


  // loop over DTOFTDCDigiHits : these are the TDC hits
  // sort them into left and right hits
  // only keep hits within the time peak
  // also keep the hodoscope planes separate 
  for (unsigned int k=0; k<TDCHits.size(); k++){
    const DTOFTDCDigiHit *hit = TDCHits[k];
    int plane = hit->plane;
    int end = hit->end;
    float time = (float)hit->time * BINTDC_2_TIME;
    TOFTDCtime->Fill(time);
    if (fabsf(time-TDCTLOC)<TDCTimeCut){
      if (end){
	TDCHitsRight[plane].push_back(hit);
      } else {
	TDCHitsLeft[plane].push_back(hit);
      }
    }
  }

  float Signum = -1.;
  vector <paddle> TOFADCPaddles[2];
  vector <SingleP> TOFADCSingles[2];
  int firsttime = 1;

  // for each hodoscope plane find matches between
  // ADC left and right data or find single hits for
  // the single ended readout paddles.
  for (int plane=0; plane<2; plane++){

    // loop over right pmts to find single ended paddle hits
    // these are paddle 22 and paddle 23
    for (unsigned int i = 0; i<ADCHitsRight[plane].size(); i++) {
      const DTOFDigiHit *hitR = ADCHitsRight[plane][i];
      int paddle = hitR->bar; 
      if ((paddle == 22) || (paddle == 23)){
	struct SingleP newsingle;
	newsingle.paddle = paddle;
	newsingle.LR = 1;
	newsingle.time = (float)hitR->pulse_time*BINADC_2_TIME ;
	newsingle.adc = (float)hitR->pulse_integral -  
	  (float)hitR->pedestal/(float)hitR->nsamples_pedestal*(float)hitR->nsamples_integral;
	newsingle.OverFlow = ADCRightOverFlow[plane][i];
	TOFADCSingles[plane].push_back(newsingle);	
      }
    }
    
    // loop over left pmts to find single ended paddle hits
    // these are paddle 22 and paddle 23    
    // for the other paddle loop over the right ones and find match
    
    for (unsigned int j = 0; j<ADCHitsLeft[plane].size(); j++) {
      const DTOFDigiHit *hit = ADCHitsLeft[plane][j];
      int paddle = hit->bar;
      if ((paddle == 22) || (paddle == 23)){ // save singles of left pmts
	struct SingleP newsingle;
	newsingle.paddle = paddle;
	newsingle.LR = 0;
	newsingle.adc = (float)hit->pulse_integral -  
	  (float)hit->pedestal/(float)hit->nsamples_pedestal*(float)hit->nsamples_integral;
	newsingle.time = (float)hit->pulse_time*BINADC_2_TIME ;
	newsingle.OverFlow = ADCLeftOverFlow[plane][j];
	TOFADCSingles[plane].push_back(newsingle);
      } else {
	
	// loop over Right adc hits find match with the left and prepare paddle hits
	for (unsigned int i = 0; i<ADCHitsRight[plane].size(); i++) {
	  const DTOFDigiHit *hitR = ADCHitsRight[plane][i];
	  if (hitR->bar == paddle){
	    struct paddle newpaddle;
	    newpaddle.paddle = paddle;
	    newpaddle.timeL = (float)hit->pulse_time*BINADC_2_TIME ;
	    newpaddle.timeR = (float)hitR->pulse_time*BINADC_2_TIME ;
	    newpaddle.mt = (newpaddle.timeL + newpaddle.timeR)/2.;
	    newpaddle.td = Signum*(newpaddle.timeL - newpaddle.timeR)/2.;
	    newpaddle.adcL = (float)hit->pulse_integral -  
	      (float)hit->pedestal/(float)hit->nsamples_pedestal*(float)hit->nsamples_integral;
	    newpaddle.adcR = (float)hitR->pulse_integral -  
	      (float)hitR->pedestal/(float)hitR->nsamples_pedestal*(float)hitR->nsamples_integral;
	    newpaddle.OverFlowL =  ADCLeftOverFlow[plane][j];
	    newpaddle.OverFlowR =  ADCRightOverFlow[plane][i];
	    if ((paddle != 22) && (paddle !=23)) {
	      //cout<<newpaddle.OverFlowL<<"   "<< newpaddle.OverFlowR <<endl;
	      TOFADCPaddles[plane].push_back(newpaddle);
	    }
	  }
	}
      }
    }
  }

  vector <paddle> TOFTDCPaddles[2];
  vector <SingleP> TOFTDCSingles[2];
  firsttime = 1;

  // now do the same thing for TDC hits
  // find matches between left and right and treat the single ended paddles separately

  for (int plane=0; plane<2; plane++){

    for (unsigned int j = 0; j<TDCHitsRight[plane].size(); j++) {
      const DTOFTDCDigiHit *hit = TDCHitsRight[plane][j];
      int paddle = hit->bar;
      if ((paddle == 22) || (paddle == 23)){
	struct SingleP newsingle;
	newsingle.paddle = paddle;
	newsingle.LR = 1;
	newsingle.time = (float)hit->time*BINTDC_2_TIME ;
	TOFTDCSingles[plane].push_back(newsingle);
      }
    }

    for (unsigned int j = 0; j<TDCHitsLeft[plane].size(); j++) {
      const DTOFTDCDigiHit *hit = TDCHitsLeft[plane][j];
      int paddle = hit->bar;
      if ((paddle == 22) || (paddle == 23)){
	struct SingleP newsingle;
	newsingle.paddle = paddle;
	newsingle.LR = 0;
	newsingle.time = (float)hit->time*BINTDC_2_TIME ;
	TOFTDCSingles[plane].push_back(newsingle);
      } else {
	for (unsigned int i = 0; i<TDCHitsRight[plane].size(); i++) {
	  const DTOFTDCDigiHit *hitR = TDCHitsRight[plane][i];
	  if (hitR->bar == paddle){
	    struct paddle newpaddle;
	    newpaddle.paddle = paddle;
	    newpaddle.timeL = (float)hit->time*BINTDC_2_TIME ;
	    newpaddle.timeR = (float)hitR->time*BINTDC_2_TIME ;
	    newpaddle.mt = (newpaddle.timeL + newpaddle.timeR)/2.;
	    newpaddle.td = Signum*(newpaddle.timeL - newpaddle.timeR)/2.; // left side get positive x
	    if ((paddle != 22) && (paddle !=23)) {
	      TOFTDCPaddles[plane].push_back(newpaddle);
	    }
	  }
	}
      }
    }
  }
 
  pthread_mutex_lock(&mutex);
  Event = eventnumber;
  TShift = TimingShift;
  Nhits = TOFTDCPaddles[0].size() + TOFTDCPaddles[1].size();
  int cnt = 0;
  int AllHits[4]={0,0,0,0};

  if (Nhits<MaxHits){
    for (unsigned int k = 0; k<TOFTDCPaddles[0].size() ; k++){
      Plane[cnt] = 0;
      Paddle[cnt] = TOFTDCPaddles[0][k].paddle;
      MeanTime[cnt] = TOFTDCPaddles[0][k].mt;
      TimeDiff[cnt] = TOFTDCPaddles[0][k].td;
      cnt++;
      AllHits[0]++;
    }
    AllHits[0] = cnt;
    for (unsigned int k = 0; k<TOFTDCPaddles[1].size() ; k++){
      Plane[cnt] = 1;
      Paddle[cnt] = TOFTDCPaddles[1][k].paddle;
      MeanTime[cnt] = TOFTDCPaddles[1][k].mt;
      TimeDiff[cnt] = TOFTDCPaddles[1][k].td;
      cnt++;
      AllHits[1]++;
    }
  } else {
    Nhits = 0;
  }

  cnt = 0;
  NhitsA = TOFADCPaddles[0].size() + TOFADCPaddles[1].size();
  if (NhitsA<MaxHits){
    for (unsigned int k = 0; k<TOFADCPaddles[0].size() ; k++){
      PlaneA[cnt] = 0;
      PaddleA[cnt] = TOFADCPaddles[0][k].paddle;
      MeanTimeA[cnt] = TOFADCPaddles[0][k].mt;
      TimeDiffA[cnt] = TOFADCPaddles[0][k].td;
      ADCL[cnt] = TOFADCPaddles[0][k].adcL;
      ADCR[cnt] = TOFADCPaddles[0][k].adcR;
      OFL[cnt] = TOFADCPaddles[0][k].OverFlowL;
      OFR[cnt] = TOFADCPaddles[0][k].OverFlowR;
      cnt++;
      AllHits[2]++;
    }

    for (unsigned int k = 0; k<TOFADCPaddles[1].size() ; k++){
      PlaneA[cnt] = 1;
      PaddleA[cnt] = TOFADCPaddles[1][k].paddle;
      MeanTimeA[cnt] = TOFADCPaddles[1][k].mt;
      TimeDiffA[cnt] = TOFADCPaddles[1][k].td;
      ADCL[cnt] = TOFADCPaddles[1][k].adcL;
      ADCR[cnt] = TOFADCPaddles[1][k].adcR;
      OFL[cnt] = TOFADCPaddles[1][k].OverFlowL;
      OFR[cnt] = TOFADCPaddles[1][k].OverFlowR;
      cnt++;
      AllHits[3]++;
    }
  } else {
    NhitsA = 0;
  }

  NsinglesA = TOFADCSingles[0].size() + TOFADCSingles[1].size();
  NsinglesT = TOFTDCSingles[0].size() + TOFTDCSingles[1].size();

  unsigned int j=0;
  for (unsigned int k = 0; k<TOFADCSingles[0].size(); k++){
    PlaneSA[k] = 0;
    PaddleSA[k] = TOFADCSingles[0][k].paddle ;
    LRA[k] = TOFADCSingles[0][k].LR ;
    ADCS[k] = TOFADCSingles[0][k].adc;
    TADCS[k] = TOFADCSingles[0][k].time; 
    OF[k] = TOFADCSingles[0][k].OverFlow;
    j++;
  }
  for (unsigned int k = 0; k<TOFADCSingles[1].size(); k++){
    PlaneSA[k+j] = 1;
    PaddleSA[k+j] = TOFADCSingles[1][k].paddle ;
    LRA[k+j] = TOFADCSingles[1][k].LR ;
    ADCS[k+j] = TOFADCSingles[1][k].adc;
    TADCS[k+j] = TOFADCSingles[1][k].time; ;
    OF[k+j] = TOFADCSingles[1][k].OverFlow;
  }
  j = 0;
  for (unsigned int k = 0; k<TOFTDCSingles[0].size(); k++){
    PlaneST[k] = 0;
    PaddleST[k] = TOFTDCSingles[0][k].paddle ;
    LRT[k] = TOFTDCSingles[0][k].LR ;
    TDCST[k] = TOFTDCSingles[0][k].time; ;
    j++;
  }
  for (unsigned int k = 0; k<TOFTDCSingles[1].size(); k++){
    PlaneST[k+j] = 1;
    PaddleST[k+j] = TOFTDCSingles[1][k].paddle ;
    LRT[k+j] = TOFTDCSingles[1][k].LR ;
    TDCST[k+j] = TOFTDCSingles[1][k].time; ;
  }

  if (((AllHits[0]>0) &&  (AllHits[1]>0)) ||
      ((AllHits[2]>0) &&  (AllHits[3]>0))){
    t3->Fill();
  }

  pthread_mutex_unlock(&mutex);

  /*
  if (!(eventnumber%10000000)){
    WriteRootFile();
  }
  */

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TOF_calib::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TOF_calib::fini(void)
{
  // Called before program exit after event processing is finished.

  WriteRootFile();

  return NOERROR;
}


jerror_t JEventProcessor_TOF_calib::WriteRootFile(void){

  sprintf(ROOTFileName,"tofdata_run%d.root",RunNumber);
  ROOTFile = new TFile(ROOTFileName,"recreate");

  pthread_mutex_lock(&mutex);
  ROOTFile->cd();
  TOFTDCtime->Write();
  TOFADCtime->Write();
  t3->Write();
  t3->AutoSave("SaveSelf");

  ROOTFile->Close();
  pthread_mutex_unlock(&mutex);
  return NOERROR;
 
}

jerror_t JEventProcessor_TOF_calib::MakeHistograms(void){

  /*
  cout<<endl;
  cout<<"CALL MakeHistograms for TOF_calib!!!! "<<endl;
  cout<<endl;
  */

  if (first){

    //cout<<"SETUP HISTOGRAMS AND TREE FOR RUN "<<RunNumber<<flush<<endl;

    first = 0;

    TOFTDCtime = new TH1F("TOFTDCtime","TOF CAEN TDC times", 8000, 0., 4000.);
    TOFADCtime = new TH1F("TOFADCtime","TOF ADC times", 800, 0., 400.);

    t3 = new TTree("t3","TOF Hits");
    t3->Branch("Event", &Event,"Event/I");

    t3->Branch("Nhits", &Nhits,"Nhits/I");
    t3->Branch("TShift",&TShift,"TShift/F");
    t3->Branch("Plane",Plane,"Plane[Nhits]/I");
    t3->Branch("Paddle",Paddle,"Paddle[Nhits]/I");
    t3->Branch("MeanTime",MeanTime,"MeanTime[Nhits]/F");
    t3->Branch("TimeDiff",TimeDiff,"TimeDiff[Nhits]/F");
  
    t3->Branch("NhitsA", &NhitsA,"NhitsA/I");
    t3->Branch("PlaneA",PlaneA,"PlaneA[NhitsA]/I");
    t3->Branch("PaddleA",PaddleA,"PaddleA[NhitsA]/I");
    t3->Branch("MeanTimeA",MeanTimeA,"MeanTimeA[NhitsA]/F");
    t3->Branch("TimeDiffA",TimeDiffA,"TimeDiffA[NhitsA]/F");
    t3->Branch("ADCL",ADCL,"ADCL[NhitsA]/F");
    t3->Branch("ADCR",ADCR,"ADCR[NhitsA]/F");
    t3->Branch("OFL",OFL,"OFL[NhitsA]/I");
    t3->Branch("OFR",OFR,"OFR[NhitsA]/I");
    
    t3->Branch("NsinglesA", &NsinglesA,"NsinglesA/I");
    t3->Branch("PlaneSA",PlaneSA,"PlaneSA[NsinglesA]/I");
    t3->Branch("PaddleSA",PaddleSA,"PaddleSA[NsinglesA]/I");
    t3->Branch("LRA",LRA,"LRA[NsinglesA]/I"); //LA=0,1 (left/right)
    t3->Branch("ADCS",ADCS,"ADCS[NsinglesA]/F");
    t3->Branch("OF",OF,"OF[NsinglesA]/I");
    t3->Branch("TADCS",TADCS,"TADCS[NsinglesA]/F");

    t3->Branch("NsinglesT", &NsinglesT,"NsinglesT/I");
    t3->Branch("PlaneST",PlaneST,"PlaneSA[NsinglesT]/I");
    t3->Branch("PaddleST",PaddleST,"PaddleSA[NsinglesT]/I");
    t3->Branch("LRT",LRT,"LRT[NsinglesT]/I"); //LA=0,1 (left/right)
    t3->Branch("TDCST",TDCST,"TDCST[NsinglesT]/F");

  }
  
  return NOERROR;
}

