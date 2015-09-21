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
    AddString(items, "module_lhit", "%d", ee.first->module);
    AddString(items, "module_rhit", "%d", ee.second->module);
    AddString(items, "t_tdc_lhit", "%f", ee.first->time_tdc);
    AddString(items, "t_tdc_rhit", "%f", ee.second->time_tdc);
    AddString(items, "t_fadc_lhit", "%f", ee.first->time_fadc);
    AddString(items, "t_fadc_rhit", "%f", ee.second->time_fadc);
    AddString(items, "integral_lhit", "%f", ee.first->integral);
    AddString(items, "integral_rhit", "%f", ee.second->integral);
    AddString(items, "pulse_peak_lhit", "%f", ee.first->pulse_peak);
    AddString(items, "pulse_peak_rhit", "%f", ee.second->pulse_peak);
  }
		
};

#endif // _DPSCPair_

