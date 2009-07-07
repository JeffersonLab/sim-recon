#ifndef _DHDDMBCALHit_
#define _DHDDMBCALHit_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;

class DHDDMBCALHit:public JObject{
	public:
		JOBJECT_PUBLIC(DHDDMBCALHit);
				
		int module;
		int layer;
		int sector;
		float E;
		float t;
		float zLocal;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "module", "%i", module);
			AddString(items, "layer", "%5d", layer);
			AddString(items, "sector", "%5d", sector);
			AddString(items, "E", "%4.3f", E);
			AddString(items, "t", "%4.3f", t);
			AddString(items, "zLocal", "%4.3f", zLocal);
		}
};

#endif // _DHDDMBCALHit_

