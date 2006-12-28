// $Id$
//
//    File: DMagneticFieldMap.h
// Created: Thu Dec 21 12:50:04 EST 2006
// Creator: davidl (on Linux alkaid 2.6.9-42.0.2.ELsmp x86_64)
//

#ifndef _DMagneticFieldMap_
#define _DMagneticFieldMap_

#include <JANA/jerror.h>

class DMagneticFieldMap{
	public:
		DMagneticFieldMap(){}
		virtual ~DMagneticFieldMap(){}
		
		virtual void GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const = 0;
		virtual void GetGradient(double x, double y, double z, double &dBdx, double &dBdy, double &dBdz, int method=0) const;

};

#endif // _DMagneticFieldMap_

