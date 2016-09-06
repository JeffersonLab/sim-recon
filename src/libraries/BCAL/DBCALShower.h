#ifndef _DBCALShower_
#define _DBCALShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <math.h>
#include <DMatrix.h>
using namespace jana;

class DBCALShower:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALShower);

    float E;
    float E_raw;
    float E_preshower;
    float x;
    float y;
    float z;
    float t;
    float xErr;
    float yErr;
    float zErr;
    float tErr;
    int N_cell;

    //for now errors are stored both in xErr,yErr,zErr and in xyzCovariance.
    //This is redundant and should be fixed.
    DMatrix ExyztCovariance;
  
	void toStrings(vector<pair<string,string> > &items)const{
	                /*Old, easier to compare r-phi rather than x-y, for Truth Hits
			AddString(items, "x", "%5.2f", x);
			AddString(items, "y", "%5.2f", y);
			*/
			AddString(items, "r", "%5.1f", sqrt(x*x+y*y));
			AddString(items, "phi", "%5.3f",atan2(y,x));
			AddString(items, "z", "%5.1f", z);
			AddString(items, "t", "%5.1f", t);
			AddString(items, "E", "%5.3f", E);
			AddString(items, "E_preshower", "%5.3f", E_preshower);
			AddString(items, "N_cell", "%d", N_cell);
	}
};

#endif // _DBCALShower_

