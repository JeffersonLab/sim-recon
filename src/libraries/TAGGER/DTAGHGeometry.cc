//
// File: DTAGHGeometry.cc
// Created: Sat Jul 5 10:18:56 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#include <stdlib.h>
#include <iostream>
#include <map>

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>
#include "DTAGHGeometry.h"

const unsigned int DTAGHGeometry::kCounterCount = 274;

//---------------------------------
// DTAGHGeometry    (Constructor)
//---------------------------------
DTAGHGeometry::DTAGHGeometry(JEventLoop *loop)
{
   /* read tagger set endpoint energy from calibdb */
   std::map<string,double> result1;
   loop->GetCalib("/PHOTON_BEAM/endpoint_energy", result1);
   if (result1.find("PHOTON_BEAM_ENDPOINT_ENERGY") == result1.end()) {
      std::cerr << "Error in DTAGHGeometry constructor: "
                << "failed to read photon beam endpoint energy "
                << "from calibdb at /PHOTON_BEAM/endpoint_energy" << std::endl;
      m_endpoint_energy_GeV = 0;
   }
   else {
      m_endpoint_energy_GeV = result1["PHOTON_BEAM_ENDPOINT_ENERGY"];
   }

   /* read hodoscope counter energy bounds from calibdb */
   std::vector<std::map<string,double> > result2;
   loop->GetCalib("/PHOTON_BEAM/hodoscope/scaled_energy_range", result2);
   if (result2.size() != kCounterCount) {
      jerr << "Error in DTAGHGeometry constructor: "
           << "failed to read photon beam scaled_energy_range table "
           << "from calibdb at /PHOTON_BEAM/hodoscope/scaled_energy_range" << std::endl;
      for (unsigned int i=0; i <= TAGH_MAX_COUNTER; ++i) {
         m_counter_xlow[i] = 0;
         m_counter_xhigh[i] = 0;
      }
   }
   else {
      for (unsigned int i=0; i < result2.size(); ++i) {
         int ctr = (result2[i])["counter"];
         m_counter_xlow[ctr] = (result2[i])["xlow"];
         m_counter_xhigh[ctr] = (result2[i])["xhigh"];
      }
   }
}

DTAGHGeometry::~DTAGHGeometry() { }

bool DTAGHGeometry::E_to_counter(double E, unsigned int &counter) const
{
   double x = E/m_endpoint_energy_GeV;
   for (counter=1; counter <= kCounterCount; ++counter) {
      if ( x >= m_counter_xlow[counter] &&
           x <= m_counter_xhigh[counter] )
      {
         return true;
      }
   }
   return false;
}

double DTAGHGeometry::getElow(unsigned int counter) const
{
   if (counter > 0 && counter <= kCounterCount)
      return m_endpoint_energy_GeV * m_counter_xlow[counter];
   else
      return 0;
}

double DTAGHGeometry::getEhigh(unsigned int counter) const
{
   if (counter > 0 && counter <= kCounterCount)
      return m_endpoint_energy_GeV * m_counter_xhigh[counter];
   else
      return 0;
}
