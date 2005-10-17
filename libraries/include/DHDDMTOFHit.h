// $Id$
//
//    File: DHDDMTOFHit.h
// Created: Mon Oct 17 15:01:51 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DHDDMTOFHit_
#define _DHDDMTOFHit_

#include "DObject.h"
#include "DFactory.h"

class DHDDMTOFHit:public DObject{
    public:
        HDCLASSDEF(DHDDMTOFHit);

        int orientation;  // 0: vertical,  1: horizontal
        int end;          // 0: left/top,  1: right/bottom
        float y;          // x/y position of bar center
        float t;          // time of light at end of bar
        float E;          // attenuated energy deposition

};

#endif // _DHDDMTOFHit_

