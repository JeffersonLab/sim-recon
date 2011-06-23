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
#include "DTOFHit.h"

/// \htmlonly
/// <A href="index.html#legend">
///	<IMG src="CORE.png" width="100">
///	</A>
/// \endhtmlonly

/// 2-plane (4-fold) TOF coincidences. The 2-hit coincidences come from DTOFHit objects
/// which are combined into coincidnces between the two planes to form 4-D space points
/// which are represented by DTOFPoint objects.

class DTOFPoint_factory:public JFactory<DTOFPoint>{
 public:
  DTOFPoint_factory(){};
  ~DTOFPoint_factory(){};

  double VELOCITY   ;
  double HALFPADDLE ;
  double BARWIDTH   ;
  
  private:
  jerror_t brun(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
  jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
 

};

#endif // _DTOFPoint_factory_

