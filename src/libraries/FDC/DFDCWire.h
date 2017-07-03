// $Id$
//
//    File: DFDCWire.h
// Created: Wed Nov 29 13:15 EST 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DFDCWire_
#define _DFDCWire_

#include <DCoordinateSystem.h>


class DFDCWire:public DCoordinateSystem{
   public:
      int layer;		///< 1-24
      int wire;		///< 1-N
      float u; // coordinate of wire in direction transverse to the wire
      float angle;	///< radians
      DVector3 angles; ///< radians
};

#endif // _DFDCWire_

