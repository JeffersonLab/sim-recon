#ifndef _DBCALShower_
#define _DBCALShower_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <math.h>
#include <DMatrix.h>
#include <TMatrixFSym.h>
using namespace jana;

class DBCALShower:public JObject{
	public:
		JOBJECT_PUBLIC(DBCALShower);

 DBCALShower():ExyztCovariance(5) {} ///< Constructor initializes matrix to 5x5

    float E;
    float E_raw;
    float E_preshower;
    float E_L2;
    float E_L3;
    float E_L4;
    float x;
    float y;
    float z;
    float t;
    float sigLong;
    float sigTrans;
    float sigTheta;
//    float sigTime;
    float rmsTime;
    int N_cell;
    int Q;
	TMatrixFSym ExyztCovariance;

    float const EErr() const { return sqrt(ExyztCovariance(0,0)); }
    float const xErr() const { return sqrt(ExyztCovariance(1,1)); }
    float const yErr() const { return sqrt(ExyztCovariance(2,2)); }
    float const zErr() const { return sqrt(ExyztCovariance(3,3)); }
	float const tErr() const { return sqrt(ExyztCovariance(4,4)); }
	float const XYcorr() const {
		if (xErr()>0 && yErr()>0) return ExyztCovariance(1,2)/xErr()/yErr();
		else return 0;
	}
	float const XZcorr() const {
		if (xErr()>0 && zErr()>0) return ExyztCovariance(1,3)/xErr()/zErr();
		else return 0;
	}
	float const YZcorr() const {
		if (yErr()>0 && zErr()>0) return ExyztCovariance(2,3)/yErr()/zErr();
		else return 0;
	}
	float const EXcorr() const {
		if (EErr()>0 && xErr()>0) return ExyztCovariance(0,1)/EErr()/xErr();
		else return 0;
	}
	float const EYcorr() const {
		if (EErr()>0 && yErr()>0) return ExyztCovariance(0,2)/EErr()/yErr();
		else return 0;
	}
	float const EZcorr() const {
		if (EErr()>0 && zErr()>0) return ExyztCovariance(0,3)/EErr()/zErr();
		else return 0;
	}
	float const XTcorr() const {
		if (xErr()>0 && tErr()>0) return ExyztCovariance(1,4)/xErr()/tErr();
		else return 0;
	}
	float const YTcorr() const {
		if (yErr()>0 && tErr()>0) return ExyztCovariance(2,4)/yErr()/tErr();
		else return 0;
	}
	float const ZTcorr() const {
		if (zErr()>0 && tErr()>0) return ExyztCovariance(3,4)/zErr()/tErr();
		else return 0;
	}
	float const ETcorr() const {
		if (EErr()>0 && tErr()>0) return ExyztCovariance(0,4)/EErr()/tErr();
		else return 0;
	}

	void toStrings(vector<pair<string,string> > &items)const{
			AddString(items, "E", "%5.3f", E);
			AddString(items, "x", "%5.2f", x);
			AddString(items, "y", "%5.2f", y);
			AddString(items, "z", "%5.1f", z);
			AddString(items, "t", "%5.1f", t);
			AddString(items, "r", "%5.1f", sqrt(x*x+y*y));
			AddString(items, "phi", "%5.3f",atan2(y,x));
			AddString(items, "E_preshower", "%5.3f", E_preshower);
                        AddString(items, "E_L2", "%5.3f", E_L2);
                        AddString(items, "E_L3", "%5.3f", E_L3);
                        AddString(items, "E_L4", "%5.3f", E_L4);
			AddString(items, "N_cell", "%d", N_cell);
			AddString(items, "Q", "%d", Q);
			AddString(items, "dE", "%5.3f", EErr());
			AddString(items, "dx", "%5.3f", xErr());
			AddString(items, "dy", "%5.3f", yErr());
			AddString(items, "dz", "%5.2f", zErr());
			AddString(items, "dt", "%5.3f", tErr());
			AddString(items, "EXcorr", "%5.3f", EXcorr());
			AddString(items, "EYcorr", "%5.3f", EYcorr());
			AddString(items, "EZcorr", "%5.3f", EZcorr());
			AddString(items, "ETcorr", "%5.3f", ETcorr());
			AddString(items, "XYcorr", "%5.3f", XYcorr());
			AddString(items, "XZcorr", "%5.3f", XZcorr());
			AddString(items, "XTcorr", "%5.3f", XTcorr());
			AddString(items, "YZcorr", "%5.3f", YZcorr());
			AddString(items, "YTcorr", "%5.3f", YTcorr());
			AddString(items, "ZTcorr", "%5.3f", ZTcorr());
			AddString(items, "sigLong", "%5.3f", sigLong);
			AddString(items, "sigTrans", "%5.3f", sigTrans);
			AddString(items, "sigTheta", "%5.3f", sigTheta);
//                        AddString(items, "sigTime", "%5.3f", sigTime);
                        AddString(items, "rmsTime", "%5.3f", rmsTime);
	}
};

#endif // _DBCALShower_

