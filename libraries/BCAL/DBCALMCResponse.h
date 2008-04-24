// $Id$
//
//    File: DBCALMCResponse.h
// Created: Thu Nov 17 09:56:05 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#ifndef _DBCALMCResponse_
#define _DBCALMCResponse_

#include "BCAL/DBCALGeometry.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DBCALMCResponse : public JObject{

public:
		JOBJECT_PUBLIC(DBCALMCResponse);
		    
    int module;
    int layer;
    int sector;
    DBCALGeometry::End end;
    float E;
    float t;
		
    int cellId;
	 
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0? "upstream":"downstream");
			AddString(items, "E(GeV)", "%d", E);
			AddString(items, "t(ns)", "%d", t);
		}
};

#endif // _DBCALMCResponse_

