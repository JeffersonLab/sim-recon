// $Id$
//
//    File: DPSCPair.h
// Created: Tue Mar 24 21:35:49 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _DPSCPair_
#define _DPSCPair_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSCHit.h"

class DPSCPair:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCPair);

  pair<const DPSCHit*,const DPSCHit*> ee;	// first:North(left); second:South(right)		
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "module_left", "%d", ee.first->module);
    AddString(items, "module_right", "%d", ee.second->module);
    AddString(items, "tl_tdc", "%f", ee.first->time_tdc);
    AddString(items, "tr_tdc", "%f", ee.second->time_tdc);
    AddString(items, "tl_fadc", "%f", ee.first->time_fadc);
    AddString(items, "tr_fadc", "%f", ee.second->time_fadc);
  }
		
};

#endif // _DPSCPair_

