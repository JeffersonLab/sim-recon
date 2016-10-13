// $Id$
//
//    File: DPSDigiHit.h
// Created: Wed Oct 15 16:46:01 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSDigiHit_
#define _DPSDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSGeometry.h"

class DPSDigiHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSDigiHit);
  
  DPSGeometry::Arm arm;   // North(left): 0, South(right): 1
  int column;
  uint32_t pulse_integral; ///< identified pulse integral as returned by FPGA algorithm
  uint32_t pulse_time;     ///< identified pulse time as returned by FPGA algorithm
  uint32_t pedestal;       ///< pedestal info used by FPGA (if any)
  uint32_t QF;             ///< Quality Factor from FPGA algorithms
  uint32_t nsamples_integral;    ///< number of samples used in integral 
  uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
  uint32_t pulse_peak;           ///<  maximum sample in pulse
		
  uint32_t datasource;           ///<  0=window raw data, 1=old(pre-Fall16) firmware, 2=Df250PulseData
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "arm", "%d", arm);
    AddString(items, "column", "%d", column);
    AddString(items, "pulse_integral", "%d", pulse_integral);
    AddString(items, "pulse_peak", "%d", pulse_peak);
    AddString(items, "pulse_time", "%d", pulse_time);
    AddString(items, "pedestal", "%d", pedestal);
    AddString(items, "QF", "%d", QF);
    AddString(items, "nsamples_integral", "%d", nsamples_integral);
    AddString(items, "nsamples_pedestal", "%d", nsamples_pedestal);
  }
};

#endif // _DPSDigiHit_

