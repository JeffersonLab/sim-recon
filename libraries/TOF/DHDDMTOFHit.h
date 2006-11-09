// $Id$
//
//    File: DHDDMTOFHit.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DHDDMTOFHit_
#define _DHDDMTOFHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DHDDMTOFHit:public JObject{

    public:
        HDCLASSDEF(DHDDMTOFHit);
	
	int plane;		// plane (0: vertical, 1: horizontal)
	int end;		// 0: north/top,  1: south/bottom
	int paddle;		// paddle number
	float t;		// time of light at end of bar
	float dE;		// attenuated energy deposition
};

#endif // _DHDDMTOFHit_

