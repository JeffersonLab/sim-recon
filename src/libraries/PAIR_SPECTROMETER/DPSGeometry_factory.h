
#ifndef _DPSGeometry_factory_
#define _DPSGeometry_factory_

#include <JANA/JFactory.h>
using namespace jana;

#include "DPSGeometry.h"

class DPSGeometry_factory:public JFactory<DPSGeometry>{
 public:
 DPSGeometry_factory(const char *tag="") :
  JFactory<DPSGeometry>(tag), factory_tag(tag) {}
  ~DPSGeometry_factory(){}

 private:
  jerror_t brun(JEventLoop *loop, int runnumber);   
  jerror_t erun(void);   
	
  std::string factory_tag;
};

#endif // _DPSGeometry_factory_
