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

class DPSCHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSCHit);
  
  int module;
  double dE;
  double t;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "module", "%4d", module);
    AddString(items, "dE (GeV)", "%f", dE);
    AddString(items, "t (ns)", "%f", t);
  }
		
};

#endif // _DPSCHit_

