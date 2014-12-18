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
      double time_tdc;
      double time_fadc;
      double npe_fadc;
      bool has_fADC,has_TDC;

      void toStrings(vector<pair<string,string> > &items)const{
        AddString(items, "counter_id", "%d", counter_id);
        AddString(items, "E(GeV)", "%f", E);
        AddString(items, "t(ns)", "%f", t);
	AddString(items, "time_tdc (ns)", "%f", time_tdc);
        AddString(items, "time_fadc(ns)", "%f", time_fadc);
        AddString(items, "npe_fadc", "%f", npe_fadc);
      }
};

#endif // _DTAGHHit_
