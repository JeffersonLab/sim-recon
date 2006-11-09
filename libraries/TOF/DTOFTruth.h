// $Id$
//
//    File: DTOFTruth.h
// Created: Mon Oct 17 13:58:02 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DTOFTruth_
#define _DTOFTruth_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFTruth:public JObject{

    public:
        HDCLASSDEF(DTOFTruth);

        int track;         //  track index
        int primary;       //  0: secondary, 1: primary
        float x, y, z;     //  true point of intersection
        float t;           //  true time

};

#endif // _DTOFTruth_

