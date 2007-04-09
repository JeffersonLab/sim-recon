// $Id$
//
//    File: DCDCWire.h
// Created: Wed Oct 18 11:39:35 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DCDCWire_
#define _DCDCWire_

#include <DCoordinateSystem.h>


class DCDCWire:public DCoordinateSystem{
	public:
		float phi;		// phi angle of wire at midplane in lab coordinates

		int ring;
		int straw;
		float stereo;
};

#endif // _DCDCWire_

