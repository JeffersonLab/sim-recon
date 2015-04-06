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

class DPSCTDCDigiHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCTDCDigiHit);
		
  // Add data members here. For example:
  int counter_id;

  uint32_t time;
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "counter_id", "%d", counter_id);
    AddString(items, "time", "%d", time);
  }
};

#endif // _DPSCTDCDigiHit_

