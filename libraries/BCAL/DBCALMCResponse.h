// $Id$
//
//    File: DBCALMCResponse.h
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALMCResponse_
#define _DBCALMCResponse_

#include "BCAL/DBCALGeometry.h"

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DBCALMCResponse : public JObject{

public:
		HDCLASSDEF(DBCALMCResponse);
		    
    int module;
    int layer;
    int sector;
    DBCALGeometry::End end;
    float E;
    float t;
		
    int cellId;
};

#endif // _DBCALMCResponse_

