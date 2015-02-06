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
  
  DPSGeometry::Arm arm;   // North: 0, South: 1
  int module;
  double dE;    // energy loss in GeV
  double t;     // hit time (best available)
  float sigma_t;  // uncertainty on t in ns
  bool has_fADC;  // true if this has an fADC hit
  bool has_TDC;   // true if this has an TDC hit


  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "arm", "%4d", arm);
    AddString(items, "module", "%4d", module);
    AddString(items, "dE(GeV)", "%f", dE);
    AddString(items, "t(ns)", "%f", t);
    AddString(items, "sigma_t", "%f", sigma_t);
    AddString(items, "has_fADC", "%d", (int)has_fADC);
    AddString(items, "has_TDC", "%d", (int)has_TDC);
  }
		
};

#endif // _DPSCHit_

