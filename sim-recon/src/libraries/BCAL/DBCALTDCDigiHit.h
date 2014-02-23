// $Id$
//
//    File: DBCALTDCDigiHit.h
// Created: Tue Aug  6 11:04:26 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALTDCDigiHit_
#define _DBCALTDCDigiHit_

#include <JANA/jerror.h>
#include <JANA/JObject.h>

#include <BCAL/DBCALGeometry.h>

class DBCALTDCDigiHit: public jana::JObject{
	public:
		JOBJECT_PUBLIC(DBCALTDCDigiHit);
	
		uint32_t module;
		uint32_t layer;
		uint32_t sector;
		DBCALGeometry::End end;
		uint32_t time;
		
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "time", "%d", time);
		}

};

#endif // _DBCALTDCDigiHit_

