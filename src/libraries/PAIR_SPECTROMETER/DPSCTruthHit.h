// $Id$

#ifndef _DPSCTruthHit_
#define _DPSCTruthHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DPSCTruthHit:public JObject{
 public:
  JOBJECT_PUBLIC(DPSCTruthHit);
  
  float dEdx;
  bool primary;
  int track;
  int itrack;
  int ptype;
  float x;
  float y;
  float z;
  float t;
  int column;
		
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "track", "%d", track);
    AddString(items, "itrack", "%d", itrack);
    AddString(items, "primary", "%d", primary);
    AddString(items, "ptype", "%d", ptype);
    AddString(items, "dEdx(MeV/cm)", "%1.3f", dEdx*1.0E3);
    AddString(items, "t", "%3.2f", t);
    AddString(items, "x", "%3.1f", x);
    AddString(items, "y", "%3.1f", y);
    AddString(items, "z", "%3.1f", z);
    AddString(items, "column", "%d", column);
  }
};

#endif // _DPSCTruthHit_

