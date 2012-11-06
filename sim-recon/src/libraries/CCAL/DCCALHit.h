// $Id$
//
//    File: DCCALHit.h
// Created: Tue Nov 30 14:54:34 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DCCALHit_
#define _DCCALHit_

#include <JANA/jerror.h>
#include <JANA/JFactory.h>

class DCCALHit:public jana::JObject{
	public:

		JOBJECT_PUBLIC(DCCALHit);

		DCCALHit(){}
		
		int row;
		int column;
		float x;
		float y;
		float E;
		float t;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "row", "%4d", row);
			AddString(items, "column", "%4d", column);
			AddString(items, "x(cm)", "%3.1f", x);
			AddString(items, "y(cm)", "%3.1f", y);
			AddString(items, "E(MeV)", "%2.3f", E*1000.0);
			AddString(items, "t(ns)", "%2.3f", t);
		}
};

#endif // _DCCALHit_

