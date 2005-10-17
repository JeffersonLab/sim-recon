// $Id$
//
//    File: DHDDMTOFTruth.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DHDDMTOFTruth_
#define _DHDDMTOFTruth_

#include "DObject.h"
#include "DFactory.h"

class DHDDMTOFTruth:public DObject{
    public:
        HDCLASSDEF(DHDDMTOFTruth);

        int orientation;   //  0: vertical, 1: horizontal
        int track;         //  track index
        int primary;       //  0: secondary, 1: primary
        float x, y, z;     //  true point of intersection
        float t;           //  true time

};

#endif // _DHDDMTOFTruth_

