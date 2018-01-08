
#ifndef _DLumi_factory_
#define _DLumi_factory_

#include <JANA/JFactory.h>
using namespace jana;

#include "DLumi.h"

class DLumi_factory:public JFactory<DLumi> {
 public:
  DLumi_factory(){}
  ~DLumi_factory(){}

 private:
  jerror_t brun(JEventLoop *loop, int32_t runnumber);
  jerror_t erun(void);
};

#endif // _DLumi_factory_
