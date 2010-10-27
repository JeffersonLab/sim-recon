// $Id$
//
//    File: DCDCHit.h
// Created: Thu Jun  9 10:22:37 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DCDCHit_
#define _DCDCHit_

#include "JANA/JObject.h"
#include "JANA/JFactory.h"

class DCDCHit:public JObject{
	public:
		JOBJECT_PUBLIC(DCDCHit);
		
		int ring;
		int straw;
		float dE;
		float t;
		float d;
		int itrack;
		int ptype;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "ring", "%d", ring);
			AddString(items, "straw", "%d", straw);
			AddString(items, "dE", "%10.4e", dE);
			AddString(items, "t", "%6.1f", t);
			AddString(items, "d", "%10.4e", d);
			AddString(items, "itrack", "%d", itrack);
			AddString(items, "ptype", "%d", ptype);
		}
};

#endif // _DCDCHit_

