// $Id$
//
//    File: DHDDMTOFHit.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//
// changes: Tue Jun 19 17:29:07 EDT 2007 B.Zihlmann
//          put north and south information to the same structure 

#ifndef _DHDDMTOFHit_
#define _DHDDMTOFHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DHDDMTOFHit:public JObject{

    public:
        HDCLASSDEF(DHDDMTOFHit);
	
	int plane;		// plane (0: vertical, 1: horizontal)
	int bar;		// bar number
        int ptype;              // GEANT particle type
	float t_north;		// time of light at end of bar
	float dE_north;		// attenuated energy deposition
	float t_south;		// time of light at end of bar
	float dE_south;		// attenuated energy deposition
	float x;                // hit location in global coordiantes
        float y;
	float z;
	float px;                // particle momentum
        float py;
	float pz;
        float E;                 // particle Energy
};

#endif // _DHDDMTOFHit_

