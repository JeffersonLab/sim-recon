// $Id$
//
//    File: DPSCHit.h
// Created: Wed Oct 15 16:45:33 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCHit_
#define _DPSCHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSGeometry.h"

class DPSCHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCHit);

  DPSGeometry::Arm arm;   // North(left): 0, South(right): 1
  int module;
  double t;
  double integral;
  double pulse_peak;
  double time_tdc;
  double time_fadc;
  double npe_fadc;
  bool has_fADC,has_TDC;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "arm", "%d", arm);
    AddString(items, "module", "%d", module);
    AddString(items, "t(ns)", "%f", t);
    AddString(items, "time_tdc(ns)", "%f", time_tdc);
    AddString(items, "time_fadc(ns)", "%f", time_fadc);
    AddString(items, "integral", "%f", integral);
    AddString(items, "pulse_peak", "%f", pulse_peak);
    AddString(items, "npe_fadc", "%f", npe_fadc);
    AddString(items, "has_fADC", "%d", (int)has_fADC);
    AddString(items, "has_TDC", "%d", (int)has_TDC);
  }
};

#endif // _DPSCHit_

