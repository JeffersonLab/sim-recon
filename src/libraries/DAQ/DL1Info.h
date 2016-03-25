// $Id$
// $HeadURL$
//
//    File: DL1Info.h


#ifndef _DL1Info_
#define _DL1Info_

#include <JANA/JObject.h>
#include <JANA/JObject.h>

using namespace jana;
using namespace std;

class DL1Info:public jana::JObject{
            public:
                JOBJECT_PUBLIC(DL1Info);
  
		uint32_t nsync;
		uint32_t trig_number;
		uint32_t live_time;
		uint32_t busy_time;
		uint32_t live_inst;
		uint32_t unix_time;
		
		vector<uint32_t> gtp_sc;
		vector<uint32_t> fp_sc;
		vector<uint32_t> gtp_rate;
		vector<uint32_t> fp_rate;
		
		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
		  AddString(items, "nsync"       , "%d" , nsync); 
		  AddString(items, "trig_number" , "%d" , trig_number); 
		  AddString(items, "live_time"   , "%d" , live_time); 
		  AddString(items, "busy_time"   , "%d" , busy_time); 
		  AddString(items, "live_inst"   , "%d" , live_inst); 
		  AddString(items, "unix_time"   , "%d" , unix_time); 
		  
		  AddString(items, "gtp_sc"    ,   "%d" ,   gtp_sc.size());
 		  AddString(items, "fp_sc"     ,   "%d" ,   fp_sc.size());
		  
		  AddString(items, "gtp_rate"   ,  "%d" ,  gtp_rate.size());	    
		  AddString(items, "fp_rate"    ,  "%d" ,  fp_rate.size());
		  
		}

		//		void toStrings(vector<pair<string,string> > &items)const{
		//			AddString(items, "rocid"      , "%d" , rocid);
		//			AddString(items, "timestamp"  , "%ld", timestamp);
		//			AddString(items, "Nmisc"      , "%d" , misc.size());
		//		}
		
};

#endif // _DL1Info_

