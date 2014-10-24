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

class DPSHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DPSHit);
		
  int column;
  double dE;
  double t;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "column", "%d", column);
    AddString(items, "dE(GeV)", "%f",dE);
    AddString(items, "t(ns)", "%f", t);
  }
};

#endif // _DPSHit_

