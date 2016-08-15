// $Id$
//
//    File: DTAGMHit_factory.h
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluey.phys.uconn.edu)
//

#ifndef _DTAGMHit_factory_
#define _DTAGMHit_factory_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DTAGMHit.h"
#include "DTAGMGeometry.h"

class DTAGMHit_factory: public jana::JFactory<DTAGMHit> {
   public:
      DTAGMHit_factory() {};
      ~DTAGMHit_factory() {};

      // config. parameter
      double DELTA_T_CLUSTER_MAX;

   private:
      jerror_t init(void);                                          ///< Called once at program start
      jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);    ///< Called everytime a new run number is detected
      jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);  ///< Called every event
      jerror_t erun(void);                                          ///< Called everytime run number changes, if brun has been called
      jerror_t fini(void);                                          ///< Called after last event of last event source has been processed

      void Reset_Data(void);
};

#endif // _DTAGMHit_factory_
