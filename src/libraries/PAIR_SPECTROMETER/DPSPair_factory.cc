// $Id$
//
//    File: DPSPair_factory.cc
// Created: Fri Mar 20 07:51:31 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "DPSPair_factory.h"
#include "DPSCPair.h"
#include "DPSHit.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPSPair_factory::init(void)
{
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSPair_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSPair_factory::evnt(JEventLoop *loop, int eventnumber)
{
  // pair spectrometer electron-positron pairs
  // require left-right PS and PSC coincidences
  vector<const DPSCPair*> cpairs;
  loop->Get(cpairs);
  if (cpairs.size()>1) {cout << "More than 1 coarse PS pair in event!!!" << endl; throw cpairs.size();}
  // First, require coincidence in PSC (trigger)
  if (cpairs.size()==0) return NOERROR;
  //
  vector<const DPSHit*> hits;
  loop->Get(hits);
  // 
  pair<const DPSHit*,const DPSHit*> ee;
  bool has_pair = false;
  double tdiff = 10.0; //ns
  if (hits.size()>1) {
    for (unsigned int i=0; i < hits.size()-1; i++) {
      for (unsigned int j=i+1; j < hits.size(); j++) {
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
  //
  DPSPair *pair = new DPSPair;
  pair->ee = ee;
  _data.push_back(pair);
  
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSPair_factory::erun(void)
{
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSPair_factory::fini(void)
{
  return NOERROR;
}

