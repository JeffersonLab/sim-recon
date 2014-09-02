//
// File: DTAGHGeometry.h
// Created: Sat Jul 5, 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#ifndef _DTAGHGeometry_
#define _DTAGHGeometry_

#include <string>

#include <JANA/JFactory.h>
#include <JANA/JObject.h>
using namespace jana;

#include "units.h"

#define TAGH_MAX_COUNTER 274

class DTAGHGeometry : public JObject {
 public:
   
   JOBJECT_PUBLIC(DTAGHGeometry);

   DTAGHGeometry(JEventLoop *loop, std::string tag, int runnumber);
   ~DTAGHGeometry();

   static const unsigned int kCounterCount;

   // counters are numbered 1..kCounterCount
   double getElow(unsigned int counter) const;
   double getEhigh(unsigned int counter) const;
   bool E_to_counter(double E, unsigned int &counter) const;

   void toStrings(vector<pair<string,string> > &items) const {
      AddString(items, "kCounterCount", "%d", kCounterCount);
   }
   
 private:
   double m_endpoint_energy_GeV;
   double m_counter_xlow[TAGH_MAX_COUNTER+1];
   double m_counter_xhigh[TAGH_MAX_COUNTER+1];
};

#endif // _DTAGHGeometry_
