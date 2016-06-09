#ifndef _Df250Scaler_
#define _Df250Scaler_

#include <JANA/JObject.h>
#include <JANA/JObject.h>

using namespace jana;
using namespace std;

class Df250Scaler:public jana::JObject{
            public:
                JOBJECT_PUBLIC(Df250Scaler);
  
                uint32_t nsync;
		uint32_t trig_number;
		uint32_t version;
		
		int crate;

		vector<uint32_t> fa250_sc;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
		  AddString(items, "nsync"       , "%d" , nsync); 
		  AddString(items, "trig_number" , "%d" , trig_number); 
		  AddString(items, "version"     , "%d" , version); 
		  AddString(items, "crate"       , "%d" , crate); 

		  AddString(items, "fa250_sc"    , "%d" , fa250_sc.size());
		}		
};

#endif // _Df250Scaler_

