// $Id$
//
//    File: DFCALHit.h
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFCALHit_
#define _DFCALHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DFCALHit:public JObject{
	
public:
	
  JOBJECT_PUBLIC(DFCALHit);
	
  DFCALHit(){}
    
  int row;
  int column;
  float x;
  float y;
  float E;
  float t;
  float intOverPeak;

  void toStrings(vector<pair<string,string> > &items) const {
    AddString(items, "row", "%4d", row);
    AddString(items, "column", "%4d", column);
    AddString(items, "x(cm)", "%3.1f", x);
    AddString(items, "y(cm)", "%3.1f", y);
    AddString(items, "E(MeV)", "%2.3f", E*1000.0);
    AddString(items, "t(ns)", "%2.3f", t);
    AddString(items, "integral over peak",  "%2.3f", intOverPeak);
  }
};

#endif // _DFCALHit_

