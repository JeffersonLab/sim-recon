#ifndef _DBCALTruthShower_
#define _DBCALTruthShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;


class DBCALTruthShower:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALTruthShower);

		int track;
		int primary;
		float phi;
		float r;
		float z;
		float t;
		float E;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "track", "%d", track);
			AddString(items, "primary", "%d", primary);
			AddString(items, "phi", "%1.3f", phi);
			AddString(items, "r", "%4.3f", r);
			AddString(items, "z", "%4.1f", z);
			AddString(items, "t", "%4.3f", t);
			AddString(items, "E", "%4.3f", E);
		}
};

#endif // _DBCALTruthShower_

