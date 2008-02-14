// $Id$
//
//    File: DBeamPhoton.h
// Created: Thu Feb 14 10:11:52 EST 2008
// Creator: davidl (on Darwin fwing-dhcp13.jlab.org 8.10.1 i386)
//

#ifndef _DBeamPhoton_
#define _DBeamPhoton_

#include <PID/DKinematicData.h>

class DBeamPhoton: public DKinematicData{
	public:
		HDCLASSDEF(DBeamPhoton);
		
		double t; ///< Time at which photon arrives at interaction vertex location
		
};



#endif // _DBeamPhoton_

