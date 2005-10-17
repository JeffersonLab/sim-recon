// $Id$
//
//    File: DTOFHit.h
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTOFHit_
#define _DTOFHit_

#include "DObject.h"
#include "DFactory.h"

class DTOFHit:public DObject{
    public:
        HDCLASSDEF(DTOFHit);

        int orientation;  // 0: vertical,  1: horizontal
        int end;          // 0: left/top,  1: right/bottom
        float y;          // x/y position of bar center
        float t;          // time of light at end of bar  (calibrated) 
        float E;          // attenuated energy deposition  (calibrated)
};

#endif // _DTOFHit_

