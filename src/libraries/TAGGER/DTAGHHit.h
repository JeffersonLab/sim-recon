// $Id$
//
// File: DTAGHHit.h
// Created: Sat Jul 5 07:49:15 EDT 2014
// Creator: jonesrt (on gluey.phys.uconn.edu)
//

#ifndef _DTAGHhit_
#define _DTAGHhit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTAGHHit:public jana::JObject{
   public:
      JOBJECT_PUBLIC(DTAGHHit);

      double E;
      double t;
      int counter_id;
      double integral;
      double pulse_peak;
      double time_tdc;
      double time_fadc;
      double npe_fadc;
      bool has_fADC, has_TDC, is_double;
      int bg = -1; //if MC, 0 for the photon that generated the event, nonzero otherwise //ignore if not MC

      void toStrings(vector<pair<string,string> > &items)const{
        AddString(items, "counter_id", "%d", counter_id);
        AddString(items, "E(GeV)", "%f", E);
        AddString(items, "t(ns)", "%f", t);
        AddString(items, "time_tdc(ns)", "%f", time_tdc);
        AddString(items, "time_fadc(ns)", "%f", time_fadc);
        AddString(items, "integral", "%f", integral);
        AddString(items, "pulse_peak", "%f", pulse_peak);
        AddString(items, "npe_fadc", "%f", npe_fadc);
        AddString(items, "has_fADC", "%d", (int)has_fADC);
        AddString(items, "has_TDC", "%d", (int)has_TDC);
        AddString(items, "is_double", "%d", (int)is_double);
        AddString(items, "bg", "%d", bg);
      }
};

#endif // _DTAGHHit_
