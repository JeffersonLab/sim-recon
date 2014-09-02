//
// File: DTAGMGeometry_factory.h
// Created: Sat Jul 5 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#ifndef _DTAGMGeometry_factory_
#define _DTAGMGeometry_factory_

#include <string>

#include "JANA/JFactory.h"
#include "DTAGMGeometry.h"

class DTAGMGeometry_factory : public JFactory<DTAGMGeometry> {
 public:
   DTAGMGeometry_factory(std::string tag="") :
      JFactory<DTAGMGeometry>(tag.c_str()), factory_tag(tag) {}
   ~DTAGMGeometry_factory(){}

 private:
   jerror_t brun(JEventLoop *loop, int runnumber);   
   jerror_t erun(void);   

   std::string factory_tag;
};

#endif // _DTAGMGeometry_factory_
