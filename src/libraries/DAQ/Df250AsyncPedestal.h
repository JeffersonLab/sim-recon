#ifndef _Df250AsyncPedestal_
#define _Df250AsyncPedestal_

#include <JANA/JObject.h>
#include <JANA/JObject.h>

using namespace jana;
using namespace std;

class Df250AsyncPedestal:public jana::JObject{
            public:
                JOBJECT_PUBLIC(Df250AsyncPedestal);
  
                uint32_t nsync;
		uint32_t trig_number;
		
		int crate;

		vector<uint32_t> fa250_ped;

		// This method is used primarily for pretty printing
		// the second argument to AddString is printf style format
		void toStrings(vector<pair<string,string> > &items)const{
		  AddString(items, "nsync"       , "%d" , nsync); 
		  AddString(items, "trig_number" , "%d" , trig_number); 
		  AddString(items, "crate"       , "%d" , crate); 

		  AddString(items, "fa250_ped"   , "%d" , fa250_ped.size());
		}		
};

#endif // _Df250AsyncPedestal_

