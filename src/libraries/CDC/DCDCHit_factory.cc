// $Id$
//
//    File: DCDCHit_factory.cc
// Created: Tue Aug  6 11:29:56 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <CDC/DCDCHit.h>

#include "DCDCHit_factory.h"

static double DIGI_THRESHOLD = -1.0e8;

//------------------
// init
//------------------
jerror_t DCDCHit_factory::init(void)
{
  gPARMS->SetDefaultParameter("CDC:DIGI_THRESHOLD",DIGI_THRESHOLD,
			      "Do not convert CDC digitized hits into DCDCHit objects"
			      " that would have q less than this");
  
  RemoveCorrelationHits = 1;
  
  gPARMS->SetDefaultParameter("CDCHit:RemoveCorrelationHits", RemoveCorrelationHits,
			      "Remove hits correlated in time with saturation hits!");
  
  RemoveCorrelationHitsCut = 1.5;
  gPARMS->SetDefaultParameter("CDCHit:RemoveCorrelationHitsCut", RemoveCorrelationHitsCut,
			      "Cut in units of 8ns bins to remove correlated hits with Saturation hits!");
  


  // Setting this flag makes it so that JANA does not delete the objects in _data.
  // This factory will manage this memory.
  SetFactoryFlag(NOT_OBJECT_OWNER);
  
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
  
  /// Read in calibration constants
  
  vector<double> cdc_timing_cuts;
  if (eventLoop->GetCalib("/CDC/timing_cut", cdc_timing_cuts)){
    LowTCut = -60.;
    HighTCut = 900.;
    jout << "Error loading /CDC/timing_cut ! set defaul values -60. and 900." << endl;
  } else {
    LowTCut = cdc_timing_cuts[0];
    HighTCut = cdc_timing_cuts[1];
    //jout<<"CDC Timing Cuts: "<<LowTCut<<" ... "<<HighTCut<<endl;
  }
  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
 
  // Clear _data vector
  _data.clear();  
  
  vector<const DCDCHit*> hits;
  loop->Get(hits, "Calib");

  if (hits.size()>4000){
    cout<<"Too many CDC hits, "<<hits.size()<< " bail!"<<endl;
    return NOERROR;
  }



  int RocID[4000];
  int Slot[4000];
  int Connector[4000];
  int Counter=0;
  double Time[4000];
  //double Amplitude[4000];
  int Max[4000];

  // loop over hits and find roc/slod/con numbers
  for (int k=0 ;k<(int) hits.size(); k++){
    const DCDCHit *hit = hits[k];
    vector <const Df125CDCPulse*> pulse;
    hit->Get(pulse);
    RocID[Counter]  = pulse[0]->rocid;
    Slot[Counter] = pulse[0]->slot;
    Connector[Counter] = pulse[0]->channel / 24;
    //Amplitude[Counter] = hit->amp;
    Time[Counter] = hit->t;
    Max[Counter] = 0;
    if (hit->QF > 1){
      Max[Counter] = 1;
    }
    Counter++;
  }
  
  int Mark4Removal[4000];
  memset(Mark4Removal, 0, 4000*sizeof(int));

  if (RemoveCorrelationHits) {
    
    for (int k=0 ;k<Counter; k++){
      
      if (Max[k]){
	
	for (int n=0 ;n<Counter; n++){
	  
	  if (n==k){
	    continue;
	  }
	  
	  if ((RocID[k] == RocID[n]) && (Slot[k] == Slot[n]) && (Connector[k] == Connector[n]) ){
	    
	    double dt = (Time[k] - Time[n])/8.; // units of samples (8ns)
	    if ( fabs(dt+3.5)<RemoveCorrelationHitsCut) {
	      Mark4Removal[n] = 1;
	      //cout<<"remove "<<hits[n]->ring<<" "<<hits[n]->straw<<endl;
	    }
	    
	  }
	  
	}
      }
    }
    
  }

  
  for (int k=0 ;k<(int) hits.size(); k++){
    const DCDCHit *hit = hits[k];

    // remove hits with small charge: not really used
    if (hit->q < DIGI_THRESHOLD) 
      continue;

    // remove hits with amplitudes 0 or less
    if ( hit->amp <= 0 ) { 
      continue;
    }

    // remove hits ouside of the timing cut
    if ( (hit->t < LowTCut) || (hit->t > HighTCut) ){
      continue;
    }

    // removed hits correclated with Saturation hit on same connector/reamp/HV-board
    if (Mark4Removal[k]){
      continue;
    }
    
    _data.push_back( const_cast<DCDCHit*>(hit) );
    
  }
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DCDCHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DCDCHit_factory::fini(void)
{
    return NOERROR;
}
