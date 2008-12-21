// $Id$
//
//    File: DVector2S.h
// Created: Fri Dec 28 07:19:32 EST 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DVector2S_
#define _DVector2S_

#include <JANA/jerror.h>

class DVector2S{
	public:
		DVector2S(){};
		DVector2S(double xx, double yy){x=xx; y=yy;}
		virtual ~DVector2S(){};
		
		inline double X(void) const {return x;}
		inline double Y(void) const {return y;}
		inline void Set(double xx, double yy){x=xx; y=yy;}
		inline double Mod(void){return sqrt(x*x + y*y);}
		inline double Phi(void){return atan2(y,x);}
		inline double Phi2pi(void){double a=atan2(y,x); return a<0.0 ? a+2*M_PI:a;}
		
		inline DVector2S& operator*=(const double &f){x*=f; y*=f; return *(this);}
		inline DVector2S& operator/=(const double &f){x/=f; y/=f; return *(this);}
		inline DVector2S& operator+=(const DVector2S &v){x+=v.X(); y+=v.Y(); return *(this);}
		inline DVector2S& operator-=(const DVector2S &v){x-=v.X(); y-=v.Y(); return *(this);}

	protected:
	
	
	private:
		double x;
		double y;
};

inline DVector2S operator*(const double &f, const DVector2S &vec){
	DVector2S s(f*vec.X(), f*vec.Y());
	return s;
}

inline DVector2S operator*(const DVector2S &vec, const double &f){
	DVector2S s(f*vec.X(), f*vec.Y());
	return s;
}

inline DVector2S operator/(const DVector2S &vec, const double &f){
	DVector2S s(vec.X()/f, vec.Y()/f);
	return s;
}

inline DVector2S operator+(const DVector2S &vec1, const DVector2S &vec2){
	DVector2S s(vec1.X()+vec2.X(), vec1.Y()+vec2.Y());
	return s;
}

inline DVector2S operator-(const DVector2S &vec1, const DVector2S &vec2){
	DVector2S s(vec1.X()-vec2.X(), vec1.Y()-vec2.Y());
	return s;
}


#endif // _DVector2S_

