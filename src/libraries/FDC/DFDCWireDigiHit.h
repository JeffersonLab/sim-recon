// $Id$
//
//    File: DFDCWireDigiHit.h
// Created: Wed Aug  7 11:54:06 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DFDCWireDigiHit_
#define _DFDCWireDigiHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

class DFDCWireDigiHit:public jana::JObject{
	public:
		JOBJECT_PUBLIC(DFDCWireDigiHit);
		
		uint32_t package;
		uint32_t chamber;
		uint32_t wire;
		uint32_t time;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "package", "%d", package);
			AddString(items, "chamber", "%d", chamber);
			AddString(items, "wire", "%d", wire);
			AddString(items, "time", "%d", time);
		}
		
};

#endif // _DFDCWireDigiHit_

