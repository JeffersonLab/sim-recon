// $Id$
//
//    File: DPSHit.h
// Created: Wed Oct 15 16:45:01 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSHit_
#define _DPSHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSGeometry.h"

class DPSHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSHit);

  DPSGeometry::Arm arm;   // North(left): 0, South(right): 1
  int column;
  double E;  
  double t;
  double integral;
  double pulse_peak;
  double npix_fadc;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "arm", "%d", arm);
    AddString(items, "column", "%d", column);
    AddString(items, "E(GeV)", "%f", E);
    AddString(items, "t(ns)", "%f", t);
    AddString(items, "integral", "%f", integral);
    AddString(items, "pulse_peak", "%f", pulse_peak);
    AddString(items, "npix_fadc", "%f", npix_fadc);
  }
};

#endif // _DPSHit_

