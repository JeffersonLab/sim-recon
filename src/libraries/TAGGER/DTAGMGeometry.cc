//
// File: DTAGMGeometry.cc
// Created: Sat Jul 5 10:18:56 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#include <stdlib.h>
#include <iostream>
#include <map>

#include <JANA/JApplication.h>
#include <JANA/JEvent.h>
#include "DTAGMGeometry.h"

const unsigned int DTAGMGeometry::kRowCount = 5;
const unsigned int DTAGMGeometry::kColumnCount = 102;
const double DTAGMGeometry::kFiberWidth = 0.2; // cm
const double DTAGMGeometry::kFiberLength = 2.0; // cm


//---------------------------------
// DTAGMGeometry    (Constructor)
//---------------------------------
DTAGMGeometry::DTAGMGeometry(JEventLoop *loop, std::string tag, int runnumber)
{
   /* read tagger set endpoint energy from calibdb */
   char dbname1[80];
   if (tag == "")
      sprintf(dbname1, "/PHOTON_BEAM/endpoint_energy:%d", runnumber);
   else
      sprintf(dbname1, "/PHOTON_BEAM/endpoint_energy:%d:%s", runnumber,
                                                         tag.c_str());
   std::map<string,double> result1;
   loop->GetCalib(dbname1, result1);
   if (result1.find("PHOTON_BEAM_ENDPOINT_ENERGY") == result1.end()) {
      std::cerr << "Error in DTAGMGeometry constructor: "
                << "failed to read photon beam endpoint energy "
                << "from calibdb at " << dbname1
                << std::endl;
      m_endpoint_energy_GeV = 0;
   }
   else {
      m_endpoint_energy_GeV = result1["PHOTON_BEAM_ENDPOINT_ENERGY"];
   }

   /* read microscope channel energy bounds from calibdb */
   char dbname2[80];
   if (tag == "") 
      sprintf(dbname2, "/PHOTON_BEAM/microscope/scaled_energy_range:%d",
                                                runnumber);
   else
      sprintf(dbname2, "/PHOTON_BEAM/microscope/scaled_energy_range:%d:%s",
                                                runnumber, tag.c_str());
   std::vector<std::map<string,double> > result2;
   loop->GetCalib(dbname2, result2);
   if (result2.size() != kColumnCount) {
      std::cerr << "Error in DTAGMGeometry constructor: "
                << "failed to read photon beam scaled_energy_range table "
                << "from calibdb at " << dbname2
                << std::endl;
      for (unsigned int i=0; i <= TAGM_MAX_COLUMN; ++i) {
         m_column_xlow[i] = 0;
         m_column_xhigh[i] = 0;
      }
   }
   else {
      for (unsigned int i=0; i < result2.size(); ++i) {
         int column = (result2[i])["column"];
         m_column_xlow[column] = (result2[i])["xlow"];
         m_column_xhigh[column] = (result2[i])["xhigh"];
      }
   }
}

DTAGMGeometry::~DTAGMGeometry() { }

bool DTAGMGeometry::E_to_column(double E, unsigned int &column) const
{
   double x = E/m_endpoint_energy_GeV;
   for (column=1; column <= kColumnCount; ++column) {
      if ( x >= m_column_xlow[column] &&
           x <= m_column_xhigh[column] )
      {
         return true;
      }
   }
   return false;
}

double DTAGMGeometry::getElow(unsigned int column)
const
{
   if (column > 0 && column <= kColumnCount)
      return m_endpoint_energy_GeV * m_column_xlow[column];
   else
      return 0;
}

double DTAGMGeometry::getEhigh(unsigned int column)
const
{
   if (column > 0 && column <= kColumnCount)
      return m_endpoint_energy_GeV * m_column_xhigh[column];
   else
      return 0;
}
