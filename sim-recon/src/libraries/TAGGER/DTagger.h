// $Id$
//
//    File: DTagger.h
// Created: Tue Aug 10 07:49:15 EDT 2010
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DTagger_
#define _DTagger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTagger:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTagger);

		double E;
		double t;
		int row;
		int column;

		void toStrings(vector<pair<string,string> > &items)const{
		  AddString(items, "row", "%d", row);
		  AddString(items, "column", "%d", column);
		  AddString(items, "E", "%f", E);
		  AddString(items, "t", "%f", t);
		}
};

#endif // _DTagger_

