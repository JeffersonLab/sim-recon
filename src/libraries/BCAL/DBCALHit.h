// $Id$
//
//    File: DBCALHit.h
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DBCALHit_
#define _DBCALHit_

#include "BCAL/DBCALGeometry.h"

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DBCALHit:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALHit);
		
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
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "E(GeV)", "%2.3f", E);
			AddString(items, "t(ns)", "%4.2f", t);
		}
};

#endif // _DBCALHit_

