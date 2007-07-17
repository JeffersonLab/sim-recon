// $Id: $
//
//    File: DMCFCALHit.h
// Created: Mon Jul 16 22:03:18 EDT 2007
// Creator: shepherd
//

#ifndef _DMCFCALHit_
#define _DMCFCALHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DMCFCALHit:public JObject{
	
public:
    
    HDCLASSDEF(DMCFCALHit);
    
    DMCFCALHit(){}
    
    int column;
    int row;
    float E;
    float t;
};

#endif // _DMCFCALHit_

