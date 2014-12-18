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
      double time_tdc;
      double time_fadc;
      double npix_fadc;
      bool has_fADC,has_TDC;

      void toStrings(vector<pair<string,string> > &items) const {
        AddString(items, "row", "%d", row);
        AddString(items, "column", "%d", column);
        AddString(items, "E(GeV)", "%f", E);
        AddString(items, "t(ns)", "%f", t);
	AddString(items, "time_tdc(ns)","%f",time_tdc);
        AddString(items, "time_fadc(ns)", "%f", time_fadc);
        AddString(items, "npix_fadc", "%f", npix_fadc);
      }
};

#endif // _DTAGMHit_
