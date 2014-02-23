#ifndef _DBCALPoint_factory_OLDSMEAR_
#define _DBCALPoint_factory_OLDSMEAR_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

using namespace jana;

#include "BCAL/DBCALPoint.h"

class DBCALHit;

class DBCALPoint_factory_OLDSMEAR : public JFactory<DBCALPoint> {

 public:

  DBCALPoint_factory_OLDSMEAR() {}
  ~DBCALPoint_factory_OLDSMEAR() {}
  const char* Tag(void){return "OLDSMEAR";}

 private:
 
  class cellHits{
  	public:
		vector<const DBCALHit*> uphits;
		vector<const DBCALHit*> dnhits;
  };
 
  jerror_t init(void);
  jerror_t evnt(JEventLoop *loop, int eventnumber);

};

#endif //_DBCALPoint_factory_OLDSMEAR_
