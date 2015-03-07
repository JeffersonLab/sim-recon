#ifndef _DTPOLHit_
#define _DTPOLHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTPOLHit:public jana::JObject{
 public:
  JOBJECT_PUBLIC(DTPOLHit);
		
  int sector;   // sector number 1-32
  int ring;     // ring number 1-24
  float dE;     // Energy loss in keV
  float t;      // best time (walk-corrected tdc)
  float t_fADC; // time from fADC
  bool has_sector; 
  bool has_ring;

  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "sector", "%d", sector);
    AddString(items, "ring", "%d", sector);
    AddString(items, "dE", "%3.3f", dE);
    AddString(items, "t", "%3.3f", t);
    AddString(items, "t_fADC", "%3.3f", t_fADC);
    AddString(items, "has_sector", "%d", (int)has_sector);
    AddString(items, "has_ring", "%d", (int)has_ring);
  }
};

#endif // _DTPOLHit_
