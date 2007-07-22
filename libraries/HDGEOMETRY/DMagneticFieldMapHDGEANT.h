// $Id$
//
//    File: DMagneticFieldMapHDGEANT.h
// Created: Thu Dec 21 14:03:31 EST 2006
// Creator: davidl (on Linux alkaid 2.6.9-42.0.2.ELsmp x86_64)
//

#ifndef _DMagneticFieldMapHDGEANT_
#define _DMagneticFieldMapHDGEANT_

#include <JANA/jerror.h>
#include "DMagneticFieldMap.h"

class DMagneticFieldMapHDGEANT:public DMagneticFieldMap{
	public:
		DMagneticFieldMapHDGEANT();
		virtual ~DMagneticFieldMapHDGEANT();

		void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method) const;

};

#endif // _DMagneticFieldMapHDGEANT_

