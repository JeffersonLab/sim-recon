// $Id$
//
//    File: DBCALTruthCell.h
// Created: Thu May  5 13:07:05 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 x86_64)
//

// This object represents the barrelEMcal->bcalCell->bcalHit
// structures from HDDM

#ifndef _DBCALTruthCell_
#define _DBCALTruthCell_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DBCALTruthCell:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DBCALTruthCell);
		
		int module;
		int layer;
		int sector;
		double E;
		double t;
		double zLocal;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%2d", module);
			AddString(items, "layer", "%2d", layer);
			AddString(items, "sector", "%1d", sector);
			AddString(items, "E", "%5.3f", E);
			AddString(items, "t", "%7.2f", t);
			AddString(items, "zLocal", "%5.1f", zLocal);
		}
		
};

#endif // _DBCALTruthCell_

