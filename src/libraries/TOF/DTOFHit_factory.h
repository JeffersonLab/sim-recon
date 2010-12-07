// $Id$
//
//    File: DTOFHit_factory.h
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTOFHit_factory_
#define _DTOFHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DTOFHit.h"

class DTOFHit_factory:public JFactory<DTOFHit>{
 public:
  DTOFHit_factory(){};
  ~DTOFHit_factory(){};
  
  double C_EFFECTIVE;
  double HALFPADDLE;

 protected:
  //jerror_t init(void);					///< Called once at program start.
  jerror_t brun(JEventLoop *eventLoop, int runnumber);	        ///< Called everytime a new run number is detected.
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  //jerror_t erun(void);					///< Called everytime run number changes, provided brun has been called.
  //jerror_t fini(void);					///< Called after last event of last event source has been processed.
};

#endif // _DTOFHit_factory_

