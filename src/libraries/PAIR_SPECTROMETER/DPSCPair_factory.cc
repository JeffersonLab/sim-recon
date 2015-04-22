// $Id$
//
//    File: DPSCPair_factory.cc
// Created: Tue Mar 24 21:35:49 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "DPSCPair_factory.h"
#include "DPSCHit.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPSCPair_factory::init(void)
{
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSCPair_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSCPair_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // coarse pair spectrometer electron-positron pairs
  // save left-right coarse PS coincidences
  vector<const DPSCHit*> hits;
  loop->Get(hits);
  //
  pair<const DPSCHit*,const DPSCHit*> ee;
  bool has_pair = false;
  double tdiff = 10.0; //ns
  if (hits.size()>1) {
    for (unsigned int i=0; i < hits.size()-1; i++) {
      for (unsigned int j=i+1; j < hits.size(); j++) {
	if (!hits[i]->has_TDC||!hits[j]->has_TDC) continue;
	if (!hits[i]->has_fADC||!hits[j]->has_fADC) continue;
	double loc_tdiff = fabs(hits[i]->t-hits[j]->t);
	if (fabs(hits[i]->arm-hits[j]->arm)==1&&loc_tdiff<tdiff) {
	  tdiff = loc_tdiff;
	  has_pair = true;
	  if (hits[i]->arm==0) {
	    ee.first = hits[i];
	    ee.second = hits[j];
	  }
	  else if (hits[i]->arm==1) {
	    ee.first = hits[j];
	    ee.second = hits[i];
	  }
	}
      }
    }
  }
  if (!has_pair) return NOERROR;
  DPSCPair *pair = new DPSCPair;
  pair->ee = ee;
  _data.push_back(pair);

  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSCPair_factory::erun(void)
{
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSCPair_factory::fini(void)
{
  return NOERROR;
}

