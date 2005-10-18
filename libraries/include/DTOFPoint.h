// $Id$
//
//    File: DTOFPoint.h
// Created: Tue Oct 18 09:50:52 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFPoint_
#define _DTOFPoint_

#include "DObject.h"
#include "DFactory.h"

class DTOFPoint:public DObject{
    public:
        HDCLASSDEF(DTOFPoint);

        float x, y, z;       //reconstructed position
        float t;             //reconstructed time
        float E;             //reconstructed energy deposition
        unsigned int nhits;  //number of hits in this point
        int hits[8];         //indices to DTOFHit (temporary form)
        float chisq;         //chisquare of this point

};

#endif // _DTOFPoint_

