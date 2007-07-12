// $Id: DTOFMCHit.h 2710 2007-06-22 15:09:48Z zihlmann $
//
//    File: DTOFMCHit.h
// Created: Mon Jul  9 16:33:20 EDT 2007
// Creator: B. Zihlmann IU
//

#ifndef _DTOFMCHit_
#define _DTOFMCHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DTOFMCHit:public JObject{
    public:
        HDCLASSDEF(DTOFMCHit);

        int orientation;  // 0: vertical,  1: horizontal
	float meantime;   
	float timediff;   
	float pos;       // hit position in paddle 
	float dpos;      // estimated uncertainty in hitposition
	float dE;        // weighted energy deposition
};

#endif // _DTOFMCHit_

