//
// File: DTAGHGeometry_factory.h
// Created: Sat Jul 5 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#ifndef _DTAGHGeometry_factory_
#define _DTAGHGeometry_factory_

#include <string>

#include "JANA/JFactory.h"
#include "DTAGHGeometry.h"

class DTAGHGeometry_factory : public JFactory<DTAGHGeometry> {
 public:
   DTAGHGeometry_factory(std::string tag="") :
      JFactory<DTAGHGeometry>(tag.c_str()), factory_tag(tag) {}
   ~DTAGHGeometry_factory(){}

 private:
   jerror_t brun(JEventLoop *loop, int runnumber);   
   jerror_t erun(void);   

   std::string factory_tag;
};

#endif // _DTAGHGeometry_factory_
