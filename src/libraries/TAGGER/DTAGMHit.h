// $Id$
//
// File: DTAGMHit.h
// Created: Sat Jul 5 07:49:15 EDT 2014
// Creator: jonesrt (on gluey.phys.uconn.edu)
//

#ifndef _DTAGMhit_
#define _DTAGMhit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTAGMHit:public jana::JObject{
   public:
      JOBJECT_PUBLIC(DTAGMHit);

      double E;
      double t;
      int row;
      int column;
      double integral;
      double pulse_peak;
      double time_tdc;
      double time_fadc;
      double npix_fadc;
      bool has_fADC,has_TDC;
      int bg = -1; //if MC, 0 for the photon that generated the event, nonzero otherwise //ignore if not MC

      void toStrings(vector<pair<string,string> > &items) const {
        AddString(items, "row", "%d", row);
        AddString(items, "column", "%d", column);
        AddString(items, "E(GeV)", "%f", (float)E);
        AddString(items, "t(ns)", "%f", (float)t);
        AddString(items, "time_tdc(ns)","%f", (float)time_tdc);
        AddString(items, "time_fadc(ns)", "%f", (float)time_fadc);
        AddString(items, "integral", "%f", (float)integral);
        AddString(items, "pulse_peak", "%f", (float)pulse_peak);
        AddString(items, "npix_fadc", "%f", (float)npix_fadc);
        AddString(items, "has_fADC", "%d", (int)has_fADC);
        AddString(items, "has_TDC", "%d", (int)has_TDC);
        AddString(items, "bg", "%d", bg);
      }
};

#endif // _DTAGMHit_
