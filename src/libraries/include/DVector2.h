// This is a replacement for TVector2 that uses sse2 instructions

#ifndef _DVector2_
#define _DVector2_

#ifndef USE_SIMD 

#include <TVector2.h>
typedef TVector2 DVector2;

#else

#include <math.h>
#include <emmintrin.h>
#include <iostream>
#include <iomanip>
#define RAD2DEG 180./M_PI
using namespace std;


class DVector2{
 public:
  DVector2(){
    vec.v=_mm_setzero_pd();
  }
  DVector2(const double a, const double b){
    vec.v=_mm_setr_pd(a,b);
  }
  DVector2(__m128d v){
    vec.v=v;
  }
  ~DVector2(){};
  
  __m128d GetV() const {return vec.v;}

  // Access by indices
  double &operator() (int row){
    return vec.d[row];
  } 
  double operator() (int row) const{
    return vec.d[row];
  }

  void Set(const double a, const double b){
    vec.v=_mm_setr_pd(a,b);
  }

  double X() const {return vec.d[0];}
  double Y() const {return vec.d[1];}

  // Vector subtraction
  DVector2 operator-(const DVector2 &v1) const{
    return DVector2(_mm_sub_pd(GetV(),v1.GetV()));
  } 
  DVector2 &operator-=(const DVector2 &v1){
    vec.v=_mm_sub_pd(GetV(),v1.GetV());
    return *this;
  }
  // Vector addition
  DVector2 operator+(const DVector2 &v1) const{
    return DVector2(_mm_add_pd(GetV(),v1.GetV()));
  }
  DVector2 &operator+=(const DVector2 &v1){
    vec.v=_mm_add_pd(GetV(),v1.GetV());
    return *this;
  }
  // division by a double 
  DVector2 &operator/=(const double c) {
    __m128d scale=_mm_set1_pd(1./c);
    vec.v=_mm_mul_pd(vec.v,scale);
    return *this;
  }
  // Dot product
  double operator*(const DVector2 &v1) const{
    return (X()*v1.X()+Y()*v1.Y());
  }
  // multiplication by a double
  DVector2 &operator*=(const double c){
    __m128d scale=_mm_set1_pd(c);
    vec.v=_mm_mul_pd(scale,GetV());
    return *this;
  }
  double Mod2() const {return (X()*X()+Y()*Y());}
  double Mod() const {return sqrt(X()*X()+Y()*Y());}
  double Phi() const {return atan2(Y(),X());}

  // Angular difference between two vectors
  double DeltaPhi(const DVector2 &v1) const { 
    double twopi=2.*M_PI;
    double dphi=Phi()-v1.Phi();
    while (dphi>=M_PI) dphi-=twopi;
    while (dphi<-M_PI) dphi+=twopi;
    return dphi;
  }
  
  // return phi angle between 0 and 2pi
  double Phi_0_2pi(double angle){
    double twopi=2.*M_PI;
    while (angle>=twopi) angle-=twopi;
    while (angle<0) angle+=twopi;
    return angle;
  }
  

  void Print(){
    cout << "DVector2 (x,y)=("<<X()<<","<<Y()<<") (rho,phi)=("<< Mod()
	 <<","<<RAD2DEG*Phi()<<")"<<endl;
  }


 private:
  union dvec{
    __m128d v;
    double d[2];
  }vec;

};

//Scaling by a double
inline DVector2 operator*(const double c,const DVector2 &v1){
  __m128d scale=_mm_set1_pd(c);
  return DVector2(_mm_mul_pd(scale,v1.GetV()));
}
inline DVector2 operator*(const DVector2 &v1,const double c){
  __m128d scale=_mm_set1_pd(c);
  return DVector2(_mm_mul_pd(scale,v1.GetV()));
}
// Division by a double 
inline DVector2 operator/(const DVector2 &v1,const double c){
  __m128d scale=_mm_set1_pd(1./c);
  return DVector2(_mm_mul_pd(v1.GetV(),scale));
}

#endif // USE_SIMD
#endif // _DVector2_
