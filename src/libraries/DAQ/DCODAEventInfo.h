// $Id$
// $HeadURL$
//
//    File: DCODAEventInfo.h
// Created: Wed Oct 15 11:35:43 EDT 2014
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _DCODAEventInfo_
#define _DCODAEventInfo_

#include <JANA/JObject.h>
using namespace jana;

class DCODAEventInfo:public JObject{
	public:
		JOBJECT_PUBLIC(DCODAEventInfo);
		
		uint32_t run_number;
		uint32_t run_type;
		uint64_t event_number;
		uint16_t event_type;
		uint64_t avg_timestamp;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "run_number"    , "%d" , run_number);
			AddString(items, "run_type"      , "%d" , run_type);
			AddString(items, "event_number"  , "%ld", event_number);
			AddString(items, "event_type"    , "%d" , event_type);
			AddString(items, "avg_timestamp" , "%ld", avg_timestamp);
		}
		
};

#endif // _DCODAEventInfo_

