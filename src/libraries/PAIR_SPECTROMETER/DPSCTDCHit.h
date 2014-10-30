// $Id$
//
//    File: DPSCTDCHit.h
// Created: Wed Oct 15 16:48:32 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCTDCHit_
#define _DPSCTDCHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSGeometry.h"

class DPSCTDCHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCTDCHit);
  
  int id;
  double t;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "id", "%4d", id);
    AddString(items, "t (ns)", "%f", t);
  }
		
};

#endif // _DPSCTDCHit_

