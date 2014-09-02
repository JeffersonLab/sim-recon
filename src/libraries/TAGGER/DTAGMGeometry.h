//
// File: DTAGMGeometry.h
// Created: Sat Jul 5, 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#ifndef _DTAGMGeometry_
#define _DTAGMGeometry_

#include <string>

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
using namespace jana;

#include "units.h"

#define TAGM_MAX_ROW     5
#define TAGM_MAX_COLUMN  102


class DTAGMGeometry : public JObject {
 public:
   
   JOBJECT_PUBLIC(DTAGMGeometry);

   DTAGMGeometry(JEventLoop *loop, std::string tag, int runnumber);
   ~DTAGMGeometry();

   static const unsigned int kRowCount;
   static const unsigned int kColumnCount;
   static const double kFiberWidth;  // cm
   static const double kFiberLength; // cm

   // columns are numbered 1..kColumnCount
   double getElow(unsigned int column) const;
   double getEhigh(unsigned int column) const;
   bool E_to_column(double E, unsigned int &column) const;

   void toStrings(vector<pair<string,string> > &items) const {
      AddString(items, "kFiberWidth", "%f cm", kFiberWidth);
      AddString(items, "kFiberLength", "%f cm", kFiberLength);
      AddString(items, "kRowCount", "%d", kRowCount);
      AddString(items, "kColumnCount", "%d", kColumnCount);
   }
   
 private:
   double m_endpoint_energy_GeV;
   double m_column_xlow[TAGM_MAX_COLUMN+1];
   double m_column_xhigh[TAGM_MAX_COLUMN+1];
};

#endif // _DTAGMGeometry_
