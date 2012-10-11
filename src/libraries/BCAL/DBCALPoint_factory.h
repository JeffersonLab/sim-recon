#ifndef _DBCALPoint_factory_
#define _DBCALPoint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"

class DBCALHit;

class DBCALPoint_factory : public JFactory<DBCALPoint> {

 public:

  DBCALPoint_factory() {}
  ~DBCALPoint_factory() {}

 private:
 
  class cellHits{
  	public:
		vector<const DBCALHit*> uphits;
		vector<const DBCALHit*> dnhits;
  };
 
  jerror_t init(void);
  jerror_t evnt(JEventLoop *loop, int eventnumber);

};

#endif //_DBCALPoint_factory_
