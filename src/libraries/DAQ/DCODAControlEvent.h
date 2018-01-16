// $Id$
// $HeadURL$
//
//    File: DCODAControlEvent.h
// Created: Sat Dec 16 07:55:18 EST 2017
// Creator: davidl (on Darwin harriet.local 13.3.0 i386)
//

#ifndef _DCODAControlEvent_
#define _DCODAControlEvent_

#include <JANA/JObject.h>
using namespace jana;

class DCODAControlEvent:public JObject{
	public:
		JOBJECT_PUBLIC(DCODAControlEvent);
		
		uint16_t event_type;
		uint32_t unix_time;
		vector<uint32_t> words;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "event_type"    , "%04x" , event_type);
			AddString(items, "unix_time"     , "%ld"  , unix_time);
			AddString(items, "Nwords"        , "%d"   , words.size());
		}
		
};

#endif // _DCODAControlEvent_

