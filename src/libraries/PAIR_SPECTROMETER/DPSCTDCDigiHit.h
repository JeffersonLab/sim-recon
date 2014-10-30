// $Id$
//
//    File: DPSCTDCDigiHit.h
// Created: Wed Oct 15 16:46:32 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCTDCDigiHit_
#define _DPSCTDCDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSGeometry.h"

class DPSCTDCDigiHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCTDCDigiHit);
		
  // Add data members here. For example:
  DPSGeometry::Arm arm;   // North: 0, South: 1
  int column;
  uint32_t time; ///< TDC time measurement
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "arm", "%d", arm==0 ? "north" : "south");
    AddString(items, "column", "%d", column);
    AddString(items, "time", "%d", time);
  }
		
};

#endif // _DPSCTDCDigiHit_

