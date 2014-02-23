// $Id$
//
//    File: DTOFTDCDigiHit.h
// Created: Wed Aug  7 09:31:00 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTOFTDCDigiHit_
#define _DTOFTDCDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DTOFTDCDigiHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DTOFTDCDigiHit);
		
		int plane;      ///< plane (0: vertical, 1: horizontal)
		int bar;        ///< bar number
		int end;        ///< left/right 0/1 or North/South 0/1
		uint32_t time;	///< hit time
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "bar", "%d", bar);
			AddString(items, "plane", "%d", plane);
			AddString(items, "end", "%d", end);
			AddString(items, "time", "%d", time);
		}
		
};

#endif // _DTOFTDCDigiHit_

