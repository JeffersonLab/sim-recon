// $Id$
//
//    File: DTrigger.h
// Created: Tue Jun  7 10:15:05 EDT 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#ifndef _DTrigger_
#define _DTrigger_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTrigger:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTrigger);
		
		bool L1fired;
		double Ebcal;
		double Efcal;
		unsigned int Nschits;
		bool SCrequired;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "L1fired", "%d", L1fired ? 1:0);
			AddString(items, "Ebcal", "%5.3f", Ebcal);
			AddString(items, "Efcal", "%5.3f", Efcal);
			AddString(items, "Nschits", "%2d", Nschits);
			AddString(items, "SCrequired", "%d", SCrequired);
		}
		
};

#endif // _DTrigger_

