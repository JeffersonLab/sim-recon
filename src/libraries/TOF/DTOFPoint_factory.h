// $Id$
//
//    File: DTOFPoint_factory.h
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFPoint_factory_
#define _DTOFPoint_factory_

#include "JANA/JFactory.h"
#include "DTOFPoint.h"
#include "DTOFHitRaw.h"

class DTOFPoint_factory:public JFactory<DTOFPoint>{
 public:
  DTOFPoint_factory(){};
  ~DTOFPoint_factory(){};
  
  private:
  jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
 

};

#endif // _DTOFPoint_factory_

