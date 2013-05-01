#ifndef _DBCALPoint_factory_NEWSMEAR_
#define _DBCALPoint_factory_NEWSMEAR_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALUnifiedHit.h"

#include <TTree.h>

class DBCALHit;

class DBCALPoint_factory_NEWSMEAR : public JFactory<DBCALPoint> {

 public:
  DBCALPoint_factory_NEWSMEAR() {}
  ~DBCALPoint_factory_NEWSMEAR() {}
  const char* Tag(void){return "NEWSMEAR";}

 private:
  class cellHits{
   public:
    vector<const DBCALUnifiedHit*> uphits;
    vector<const DBCALUnifiedHit*> dnhits;
  };
 
  jerror_t evnt(JEventLoop *loop, int eventnumber);
};

#endif //_DBCALPoint_factory_NEWSMEAR_
