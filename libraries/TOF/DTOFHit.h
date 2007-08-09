// $Id$
//
//    File: DTOFHit.h
// Created: Thu Jun  9 10:05:21 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTOFHit_
#define _DTOFHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFHit:public JObject{
    public:
        HDCLASSDEF(DTOFHit);

        int orientation;  // 0: vertical,  1: horizontal
	int bar;          // bar number
        float t_north;    // time of light at end of bar  (calibrated) 
        float E_north;    // attenuated energy deposition  (calibrated)
        float t_south;    // time of light at end of bar  (calibrated) 
        float E_south;    // attenuated energy deposition  (calibrated)

        float meantime;   // equivalent to time of flight
        float timediff;    // north - south time difference
        float pos;        // hit position in paddle
        float dpos;       // estimated uncertainty in hitposition
        float dE;         // weighted energy deposition

};

#endif // _DTOFHit_

