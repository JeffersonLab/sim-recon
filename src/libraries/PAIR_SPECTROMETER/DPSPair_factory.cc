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
#include "DPSHit.h"
using namespace jana;

inline bool DPSPair_SortByTimeDifference(const DPSPair* pair1, const DPSPair* pair2)
{
  double tdiff1 = fabs(pair1->ee.first->t-pair1->ee.second->t);
  double tdiff2 = fabs(pair2->ee.first->t-pair2->ee.second->t);
  return (tdiff1<tdiff2);
}

//------------------
// init
//------------------
jerror_t DPSPair_factory::init(void)
{
  DELTA_T_PAIR_MAX = 10.0; // ns
  gPARMS->SetDefaultParameter("PSPair:DELTA_T_PAIR_MAX",DELTA_T_PAIR_MAX,
			      "Maximum difference in ns between a pair of hits"
			      " in left and right arm of fine PS");
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
  // get fine pair spectrometer hits
  vector<const DPSHit*> hits;
  loop->Get(hits);
  // form PS left-right hit pairs and sort by time difference
  pair<const DPSHit*,const DPSHit*> ee;
  if (hits.size()>1) {
    for (unsigned int i=0; i < hits.size()-1; i++) {
      for (unsigned int j=i+1; j < hits.size(); j++) {
	if (fabs(hits[i]->arm-hits[j]->arm)==1&&fabs(hits[i]->t-hits[j]->t)<DELTA_T_PAIR_MAX) {
	  if (hits[i]->arm==0) {
	    ee.first = hits[i];
	    ee.second = hits[j];
	  }
	  else if (hits[i]->arm==1) {
	    ee.first = hits[j];
	    ee.second = hits[i];
	  }
	  DPSPair *pair = new DPSPair;
	  pair->ee = ee;
	  _data.push_back(pair);
	}
      }
    }
  }
  sort(_data.begin(),_data.end(),DPSPair_SortByTimeDifference);

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

