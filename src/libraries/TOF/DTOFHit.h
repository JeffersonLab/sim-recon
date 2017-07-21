// $Id: DTOFHit.h Tue Jan 18 16:15:26 EST 2011
//
/// File:    DTOFHit.h
/// Created: Tue Jan 18 16:15:26 EST 2011
/// Creator: B. Zihlmann
/// Purpose: Container class to hold Monte Carlo data, unsmeared and 
///          smeared with the MC tag.
//

#ifndef _DTOFHit_
#define _DTOFHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFHit:public jana::JObject{
  
 public:
  JOBJECT_PUBLIC(DTOFHit);
  
  int plane;      // plane (0: vertical, 1: horizontal)
  int bar;        // bar number
  int end;        // left/right 0/1 or North/South 0/1
  float dE;       // attenuated energy deposition
  float Amp;      // Signal max Amplitude, ADC counts
  float t_fADC;   // time from adc
  float t_TDC;  // time from tdc
  float t; // walk corrected time
  bool has_fADC;
  bool has_TDC;
  
  void toStrings(vector<pair<string,string> > &items)const{
    AddString(items, "bar", "%d", bar);
    AddString(items, "plane", "%d", plane);
    AddString(items, "end", "%d", end);
    AddString(items, "dE", "%12.4f", dE);
    AddString(items, "Amp", "%12.4f", Amp);
    AddString(items, "t", "%12.4f", t);
    AddString(items, "t_TDC","%12.4f",t_TDC);
    AddString(items, "t_fADC","%12.4f",t_fADC);
  }
};

#endif // _DTOFHit_

