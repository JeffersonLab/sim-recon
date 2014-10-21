#ifndef _DBCALTruthShower_
#define _DBCALTruthShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
using namespace jana;


class DBCALTruthShower:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALTruthShower);

		int track; ///< This is the unique number that GEANT has assigned the particle
		int itrack; ///< This is the index within the MCThrown structure of this track
		int ptype; ///< This is the particle ID number
		int primary;
		float phi;
		float r;
		float z;
		float t;
		float E;
		float px;
		float py;
		float pz;

		void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "ptype", "%d", ptype);
			AddString(items, "track", "%d", track);
			AddString(items, "itrack", "%d", itrack);
			AddString(items, "primary", "%d", primary);
			AddString(items, "phi", "%1.3f", phi);
			AddString(items, "r", "%4.3f", r);
			AddString(items, "z", "%4.1f", z);
			AddString(items, "t", "%4.3f", t);
			AddString(items, "p", "%4.3f", sqrt(px*px + py*py + pz*pz));
			AddString(items, "E", "%4.3f", E);
		}
};

#endif // _DBCALTruthShower_

