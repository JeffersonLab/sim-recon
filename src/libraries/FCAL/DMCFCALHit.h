// $Id: $
//
//    File: DMCFCALHit.h
// Created: Mon Jul 16 22:03:18 EDT 2007
// Creator: shepherd
//

#ifndef _DMCFCALHit_
#define _DMCFCALHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DMCFCALHit:public JObject{
	
public:
    
    JOBJECT_PUBLIC(DMCFCALHit);
    
    DMCFCALHit(){}
    
    int column;
    int row;
    float E;
    float t;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "column", "%d", column);
			AddString(items, "row", "%d", row);
			AddString(items, "E(GeV)", "%3.2f", E);
			AddString(items, "t(ns)", "%3.2f", t);
		}
};

#endif // _DMCFCALHit_

