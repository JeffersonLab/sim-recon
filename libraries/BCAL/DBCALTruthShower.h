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
			AddString(items, "x", "%d", track);
			AddString(items, "y", "%4.3f", phi);
			AddString(items, "z", "%3.1f", r);
			AddString(items, "t", "%4.1f", z);
			AddString(items, "E", "%4.3f", E);
			AddString(items, "E", "%4.3f", t);
		}
};

#endif // _DBCALTruthShower_

