// $Id$
//
//    File: DSCTDCDigiHit.h
// Created: Tue Aug  6 13:02:22 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DSCTDCDigiHit_
#define _DSCTDCDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DSCTDCDigiHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DSCTDCDigiHit);
		
		// Add data members here. For example:
		int sector;		///< sector number 1-24
		
		uint32_t time;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "sector", "%d", sector);
			AddString(items, "time", "%d", time);
		}
		
};

#endif // _DSCTDCDigiHit_

