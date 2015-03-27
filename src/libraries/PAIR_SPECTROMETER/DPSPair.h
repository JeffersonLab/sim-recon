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
    AddString(items, "column_left", "%d", ee.first->column);
    AddString(items, "column_right", "%d", ee.second->column);
    AddString(items, "Epair", "%f", ee.first->E+ee.second->E);
    AddString(items, "El", "%f", ee.first->E);
    AddString(items, "Er", "%f", ee.second->E);
    AddString(items, "tl", "%f", ee.first->t);
    AddString(items, "tr", "%f", ee.second->t);
  }
		
};

#endif // _DPSPair_

