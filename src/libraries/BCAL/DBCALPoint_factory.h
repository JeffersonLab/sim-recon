#ifndef _DBCALPoint_factory_
#define _DBCALPoint_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"

class DBCALPoint_factory : public JFactory<DBCALPoint> {

 public:

  DBCALPoint_factory() {}
  ~DBCALPoint_factory() {}

 private:

  jerror_t evnt(JEventLoop *loop, int eventnumber);

};

#endif //_DBCALPoint_factory_
