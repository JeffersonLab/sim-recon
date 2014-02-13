// $Id$
//
//    File: DBCALTDCHit.h
// Created: Thu Aug  2 11:40:16 EDT 2012
// Creator: davidl (on Darwin eleanor.jlab.org 10.8.0 i386)
//

#ifndef _DBCALTDCHit_
#define _DBCALTDCHit_

#include <JANA/jerror.h>
#include <BCAL/DBCALGeometry.h>

class DBCALTDCHit:public JObject{

	/// This class holds data originating from the F1TDC
	/// modules connected to the BCAL

	public:
		JOBJECT_PUBLIC(DBCALTDCHit);
		
		int module;
		int layer;
		int sector;
		DBCALGeometry::End end;
		float t;

		int cellId;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%d", module);
			AddString(items, "layer", "%d", layer);
			AddString(items, "sector", "%d", sector);
			AddString(items, "end", "%s", end==0 ? "upstream":"downstream" );
			AddString(items, "t(ns)", "%4.2f", t);
		}
};

#endif // _DBCALTDCHit_

