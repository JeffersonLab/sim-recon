#ifndef _DBCALShower_
#define _DBCALShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <math.h>
using namespace jana;

class DBCALShower:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALShower);

    float E;
    float E_raw;
    float x;
    float y;
    float z;
    float t;
    float xErr;
    float yErr;
    float zErr;
    float tErr;
    int N_cell;
  
	void toStrings(vector<pair<string,string> > &items)const{
	                /*Old, easier to compare r-phi rather than x-y, for Truth Hits
			AddString(items, "x", "%5.2f", x);
			AddString(items, "y", "%5.2f", y);
			*/
	      AddString(items, "r", "%5.2f", sqrt(x*x+y*y));
			AddString(items, "phi", "%5.2f",atan2(y,x));
			AddString(items, "z", "%5.2f", z);
			AddString(items, "t", "%5.2f", t);
			AddString(items, "E", "%5.2f", E);
			AddString(items, "N_cell", "%d", N_cell);
	}
};

#endif // _DBCALShower_

