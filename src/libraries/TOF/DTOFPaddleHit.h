// $Id$
//
//    File: DTOFPaddleHit.h
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTOFPaddleHit_
#define _DTOFPaddleHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DTOFPaddleHit:public JObject{
 public:
  JOBJECT_PUBLIC(DTOFPaddleHit);
  
  int orientation;  // 0: vertical,  1: horizontal
  int bar;          // bar number
  float t_north;    // time of light at end of bar  (calibrated) 
  float E_north;    // attenuated energy deposition  (calibrated)
  float t_south;    // time of light at end of bar  (calibrated) 
  float E_south;    // attenuated energy deposition  (calibrated)
  
  float meantime;   // equivalent to time of flight
  float timediff;    // north - south time difference
  float pos;        // hit position in paddle
  float dpos;       // estimated uncertainty in hitposition
  float dE;         // weighted energy deposition
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "orientation", "%d", orientation);
    AddString(items, "bar", "%d", bar);
    AddString(items, "t_north", "%1.3f", t_north);
    AddString(items, "E_north", "%1.3f", E_north);
    AddString(items, "t_south", "%1.3f", t_south);
    AddString(items, "E_south", "%1.3f", E_south);
    AddString(items, "meantime", "%1.3f", meantime);
    AddString(items, "timediff", "%1.3f", timediff);
    AddString(items, "pos", "%1.3f", pos);
    AddString(items, "dpos", "%1.3f", dpos);
    AddString(items, "dE", "%1.3f", dE);
  }
};

#endif // _DTOFPaddleHit_

