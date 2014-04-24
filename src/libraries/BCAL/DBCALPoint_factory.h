#ifndef _DBCALPoint_factory_
#define _DBCALPoint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"

#include <TTree.h>

class DBCALHit;

class DBCALPoint_factory : public JFactory<DBCALPoint> {

 public:
  DBCALPoint_factory() {}
  ~DBCALPoint_factory() {}

 private:
  class cellHits{
   public:
    vector<const DBCALUnifiedHit*> uphits;
    vector<const DBCALUnifiedHit*> dnhits;
  };

  double m_z_target_center;
 
  jerror_t brun(JEventLoop *loop, int runnumber);
  jerror_t evnt(JEventLoop *loop, int eventnumber);
};

#endif //_DBCALPoint_factory_
