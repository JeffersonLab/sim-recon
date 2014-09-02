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
      double time_fadc;
      double npix_fadc;

      void toStrings(vector<pair<string,string> > &items) const {
        AddString(items, "row", "%d", row);
        AddString(items, "column", "%d", column);
        AddString(items, "E(GeV)", "%f", E);
        AddString(items, "t(ns)", "%f", t);
        AddString(items, "time_fadc(ns)", "%f", E);
        AddString(items, "npix_fadc", "%f", t);
      }
};

#endif // _DTAGMHit_
