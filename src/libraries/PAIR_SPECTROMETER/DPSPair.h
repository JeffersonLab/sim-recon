// $Id$
//
//    File: DPSPair.h
// Created: Fri Mar 20 07:51:31 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _DPSPair_
#define _DPSPair_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include "DPSHit.h"

class DPSPair:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSPair);
  
  pair<const DPSHit*,const DPSHit*> ee;	// first:North(left); second:South(right)	

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "column_lhit", "%d", ee.first->column);
    AddString(items, "column_rhit", "%d", ee.second->column);
    AddString(items, "E_pair", "%f", ee.first->E+ee.second->E);
    AddString(items, "E_lhit", "%f", ee.first->E);
    AddString(items, "E_rhit", "%f", ee.second->E);
    AddString(items, "t_lhit", "%f", ee.first->t);
    AddString(items, "t_rhit", "%f", ee.second->t);
    AddString(items, "integral_lhit", "%f", ee.first->integral);
    AddString(items, "integral_rhit", "%f", ee.second->integral);
    AddString(items, "pulse_peak_lhit", "%f", ee.first->pulse_peak);
    AddString(items, "pulse_peak_rhit", "%f", ee.second->pulse_peak);
  }
		
};

#endif // _DPSPair_

