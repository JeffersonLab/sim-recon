// $Id$
//
//    File: DTAGMDigiHit.h
// Created: Tue Aug  2 12:23:55 EDT 2014
// Creator: jonesrt (on Linux gluey.phys.uconn.edu)
//

#ifndef _DTAGMDigiHit_
#define _DTAGMDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTAGMDigiHit:public jana::JObject{
   public:
      JOBJECT_PUBLIC(DTAGMDigiHit);
      
      int row;
      int column;
      uint32_t pulse_integral; ///< identified pulse integral as returned by FPGA algorithm
      uint32_t pulse_time;     ///< identified pulse time as returned by FPGA algorithm
      uint32_t pedestal;       ///< pedestal info used by FPGA (if any)
      uint32_t QF;             ///< Quality Factor from FPGA algorithms
      uint32_t nsamples_integral;    ///< number of samples used in integral 
      uint32_t nsamples_pedestal;    ///< number of samples used in pedestal
      
      // This method is used primarily for pretty printing
      // the second argument to AddString is printf style format
      void toStrings(vector<pair<string,string> > &items)const{
         AddString(items, "row", "%4d", row);
         AddString(items, "column", "%4d", column);
         AddString(items, "pulse_integral", "%d", pulse_integral);
         AddString(items, "pulse_time", "%d", pulse_time);
         AddString(items, "pedestal", "%d", pedestal);
         AddString(items, "QF", "%d", QF);
      }
      
};

#endif // _DTAGMDigiHit_

