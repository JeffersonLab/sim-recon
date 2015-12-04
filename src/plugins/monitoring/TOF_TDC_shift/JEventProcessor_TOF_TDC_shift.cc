#include <stdint.h>
#include <iomanip>
#include <vector>
#include <TMath.h>
#include "JEventProcessor_TOF_TDC_shift.h"
#include <JANA/JApplication.h>
using namespace std;
using namespace jana;

#include "TOF/DTOFHit.h"
#include "TOF/DTOFDigiHit.h"
#include "TOF/DTOFTDCDigiHit.h"
#include "TTAB/DTranslationTable.h"

#include "DAQ/DCODAROCInfo.h"

// Define some constants
const float_t   ADC_BIN = 0.0625; // fADC250 pulse time resolution (ns)
const float_t   TDC_BIN = 0.0234375; // CAEN TDC time resolution (ns)

//----------------------------------------------------------------------------------
// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TOF_TDC_shift());
  }
}
//----------------------------------------------------------------------------------
JEventProcessor_TOF_TDC_shift::JEventProcessor_TOF_TDC_shift() {
}
//----------------------------------------------------------------------------------
JEventProcessor_TOF_TDC_shift::~JEventProcessor_TOF_TDC_shift() {
}
//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TOF_TDC_shift::init(void) {
  
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  // Create root folder for ST and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("TOF_TDC_shift")->cd();

  // TI remainder vs (ADC time - TDC time)
  hrocTimeRemainder_AdcTdcTimeDiff = new TH2I("hrocTimeRemainder_AdcTdcTimeDiff",";t_{ADC} - t_{TDC};TI % 6",4000,-1500,500,6,-0.5,5.5);

  // cd back to main directory
  main->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}
//----------------------------------------------------------------------------------
jerror_t JEventProcessor_TOF_TDC_shift::brun(JEventLoop *eventLoop, int32_t runnumber) {

  char filename[200];
  sprintf(filename,"TOF_TDC_shift_%6.6d.txt",runnumber);
  OUTPUT.open(filename);
  // OUTPUT << setw(6) << runnumber;

  // This is called whenever the run number changes
  return NOERROR;
}
//----------------------------------------------------------------------------------
jerror_t JEventProcessor_TOF_TDC_shift::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {

  // Get all data objects first so we minimize the time we hold the ROOT mutex lock

  // Each detector's hits
  vector<const DTOFDigiHit*>        dtofdigihits;        // TOF DigiHits
  vector<const DTOFTDCDigiHit*>     dtoftdcdigihits;     // TOF TDC DigiHits
  vector<const DCODAROCInfo*>       dcodarocinfo;        // DCODAROCInfo

  // TOF
  eventLoop->Get(dtofdigihits);
  eventLoop->Get(dtoftdcdigihits);
  eventLoop->Get(dcodarocinfo);

  Int_t TriggerBIT = -1;
  for(UInt_t i=0;i<dcodarocinfo.size();i++){
    Int_t rocid = dcodarocinfo[i]->rocid;
    // We can use any roc ID besides 1,
    // just use 11.
    if(rocid == 11){
      ULong64_t rocTime = dcodarocinfo[i]->timestamp;
      TriggerBIT = rocTime % 6;
      break;
    }
  }

  // Lock ROOT mutex so other threads won't interfere 
  japp->RootWriteLock();


  // Fill histogram of TI % 6 vs (ADC time - TDC time)
  for(UInt_t tof = 0;tof<dtofdigihits.size();tof++){
    for(UInt_t tof_TDC = 0;tof_TDC<dtoftdcdigihits.size();tof_TDC++){
      // Don't bother with matching of ADC and TDC hits,
      // just fill for all combinations

      // Make sure TI remainder is there
      if(TriggerBIT != -1){
	Double_t adc_time = dtofdigihits[tof]->pulse_time * ADC_BIN;
	Double_t tdc_time = (Double_t)dtoftdcdigihits[tof_TDC]->time * TDC_BIN;
	Double_t diff     = adc_time - tdc_time;
	// cout << adc_time << "   " << tdc_time << "   " << diff << endl;
	hrocTimeRemainder_AdcTdcTimeDiff->Fill(diff, TriggerBIT);
      }
    }
  }
  // Lock ROOT mutex so other threads won't interfere 
  japp->RootUnLock();

  return NOERROR;
}
//----------------------------------------------------------------------------------
jerror_t JEventProcessor_TOF_TDC_shift::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}
//----------------------------------------------------------------------------------
jerror_t JEventProcessor_TOF_TDC_shift::fini(void) {

  // calculate mean timing for 6 possible values of
  // TI remainder. The TI remainder bin with the highest
  // mean gives the shift value for this run.
  Double_t mean[6];
  TH1I *hproj;
  char hname[200];
  Double_t min = +99999;
  Int_t shift = -1;
  gDirectory->cd("TOF_TDC_shift");
  for(Int_t i=0;i<6;i++){
    // Get projection of (ADC time - TDC time) for each value of TI remainder
    sprintf(hname,"TOF_TDC_shift/h%d",i);
    hproj = (TH1I*)hrocTimeRemainder_AdcTdcTimeDiff->ProjectionX(hname,i+1,i+1);
    mean[i] = hproj->GetMean();
    cout << "TI remainder = " << i << " mean = " << mean[i] << endl;
    if(mean[i] < min){
      min = mean[i];
      shift = i;
    }
  }

  // shift value can only be 1, 3, or 5.
  if(!(shift == 1 || shift == 3 || shift == 5)) shift = -1;

  OUTPUT << shift << endl;

  // Called before program exit after event processing is finished.
  return NOERROR;
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
