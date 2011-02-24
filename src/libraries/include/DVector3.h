// $Id$
//
//    File: DVector3.h
//
// This header file is a replacement of TVector3 using sse2 SIMD instructions
//

#ifndef _DVector3_
#define _DVector3_

#ifndef USE_SSE2

#include <TVector3.h>
typedef TVector3 DVector3;

#else
#include <math.h>
#include <emmintrin.h>
#include <iostream>
#include <iomanip>
#define RAD2DEG 180./M_PI
using namespace std;

class DVector3{
 public:
  DVector3(){
    xy.v=_mm_setr_pd(0.,0.);
    zx.v=_mm_setr_pd(0.,0.);
  };
  DVector3(const double x,const double y,const double z){
    xy.v=_mm_setr_pd(x,y);
    zx.v=_mm_setr_pd(z,x); 
  };
  DVector3(__m128d v1,__m128d v2){
    xy.v=v1;
    zx.v=v2;
  }
  ~DVector3(){};

  // Routines to set the components of the vector
  void SetXY(const double x,const double y){
    xy.v=_mm_setr_pd(x,y);
  }
  void SetXYZ(const double x,const double y,const double z){
    xy.v=_mm_setr_pd(x,y);
    zx.v=_mm_setr_pd(z,x);
  }
  void SetMagThetaPhi(const double p, const double theta, const double phi){
    double my_p=fabs(p);
    double pt=my_p*sin(theta); 
    xy.d[0]=zx.d[1]=pt*cos(phi);
    xy.d[1]=pt*sin(phi);
    zx.d[0]=my_p*cos(theta);
  }
  // Set phi keeping theta and magnitude fixed
  void SetPhi(const double phi){
    double pt=Perp();
    xy.d[0]=zx.d[1]=pt*cos(phi);
    xy.d[1]=pt*sin(phi);
  }
  void SetX(const double x){
    xy.d[0]=zx.d[1]=x;
  } 
  void SetY(const double y){
    xy.d[1]=y;
  } 
  void SetZ(const double z){
    zx.d[0]=z;
  }

  // Get component by index
  double operator () (int ind) const {
    switch (ind){
    case 0:
      return x();
      break;
    case 1:
      return y();
      break;
    case 2:
      return z();
    default:
      // cerr << "Invalid index." <<endl;
      break;
    }
    return 0.;
  }
  
  const double x() const {return xy.d[0];};
  const double y() const {return xy.d[1];};
  const double z() const {return zx.d[0];}; 
  const double X() const {return xy.d[0];};
  const double Y() const {return xy.d[1];};
  const double Z() const {return zx.d[0];};
  const double Px() const {return xy.d[0];};
  const double Py() const {return xy.d[1];};
  const double Pz() const {return zx.d[0];}; 

  double CosTheta() const{
    double r=Mag();
    return r == 0.0 ? 1.0 : z()/r;
  }
  double Theta() const{
    return acos(CosTheta());
  }
  double Phi() const{
    return atan2(y(),x());
  }

  double Perp2() const{
    return (xy.d[0]*xy.d[0]+xy.d[1]*xy.d[1]);
  }

  double Perp() const{
    return (sqrt(Perp2()));
  }
  double Pt() const{ return Perp();};

  double Mag2() const{
    return(xy.d[0]*xy.d[0]+xy.d[1]*xy.d[1]+zx.d[0]*zx.d[0]);
  }

  double Mag() const{
    return (sqrt(Mag2()));
  };

  // Cross product of "this"=(x,y,z) and v1
  //  |x'| |y*vz-z*vy|
  //  |y'|=|z*vx-x*vz|
  //  |z'| |x*vy-y*vx|
  DVector3 Cross(const DVector3 &v1) const {
    __m128d yz=_mm_setr_pd(y(),z());
    __m128d yz1=_mm_setr_pd(v1.y(),v1.z());
    return DVector3(_mm_sub_pd(_mm_mul_pd(yz,v1.zx.v),
			       _mm_mul_pd(zx.v,yz1)),
		    _mm_sub_pd(_mm_mul_pd(xy.v,yz1),
                               _mm_mul_pd(yz,v1.xy.v))  
		    );
  }

  // Create a vector orthogonal to "this" 
  DVector3 Orthogonal() const{  
    double xx= x()<0.0 ? -x() : x();
    double yy= y()<0.0 ? -y() : y();
    double zz= z()<0.0 ? -z() : z();
    if (xx < yy) {
      return xx < zz ? DVector3(0.,z(),-y()) : DVector3(y(),-x(),0.);
    } else {
      return yy < zz ? DVector3(-z(),0.,x()) : DVector3(y(),-x(),0.);
    }
  };

  // Set the magnitude of the vector to c
  DVector3 SetMag(const double c){
    __m128d scale=_mm_set1_pd(c/Mag());
    xy.v=_mm_mul_pd(xy.v,scale);
    zx.v=_mm_mul_pd(zx.v,scale);
    return *this;
  }

  // Check for equality
  bool operator==(const DVector3 &v1) const {
    return ((x()==v1.x() && y()==v1.y() && z()==v1.z())? true : false);
  }
  // Check for inequality 
  bool operator!=(const DVector3 &v1) const{
    return ((x()!=v1.x() || y()!=v1.y() || z()!=v1.z())? true : false);
  }


  // Assignment operator
  DVector3 &operator=(const DVector3 &v1){
    xy.v=v1.xy.v; 
    zx.v=v1.zx.v;
    return *this;
  };

  // Addition
  DVector3 &operator+=(const DVector3 &v1){
    xy.v=_mm_add_pd(xy.v,v1.xy.v);
    zx.v=_mm_add_pd(zx.v,v1.zx.v);
    return *this;
  };
 
  // Subtraction
  DVector3 &operator-=(const DVector3 &v1){
    xy.v=_mm_sub_pd(xy.v,v1.xy.v);
    zx.v=_mm_sub_pd(zx.v,v1.zx.v);
    return *this;
  };

  // Add two vectors
  DVector3 operator+(const DVector3 &v2) const{
    return DVector3(_mm_add_pd(GetVxy(),v2.GetVxy()),
		    _mm_add_pd(GetVzx(),v2.GetVzx()));
  };
  // Subtract two vectors
  DVector3 operator-(const DVector3 &v2) const{
    return DVector3(_mm_sub_pd(GetVxy(),v2.GetVxy()),
		    _mm_sub_pd(GetVzx(),v2.GetVzx()));
  };



  // Scaling 
  DVector3 &operator*=(const double c){  
    __m128d scale=_mm_set1_pd(c);
    xy.v=_mm_mul_pd(GetVxy(),scale);
    zx.v=_mm_mul_pd(GetVzx(),scale);
    return *this;
  }
  
  // Unary minus
  DVector3 operator-() const{
    __m128d zero=_mm_set1_pd(0.);
    return DVector3(_mm_sub_pd(zero,xy.v),_mm_sub_pd(zero,zx.v));
  }

  // Rotate by angle a about the z-axis
  DVector3 RotateZ(const double a){
    __m128d sa=_mm_set1_pd(sin(a));
    __m128d ca=_mm_set1_pd(cos(a));
    xy.v=_mm_add_pd(_mm_mul_pd(ca,xy.v),
		    _mm_mul_pd(sa,_mm_setr_pd(-y(),x())));
    zx.d[1]=x();
    return *this;
  } 

  // Rotate by angle a about the x-axis
  DVector3 RotateX(const double a){
    __m128d sa=_mm_set1_pd(sin(a));
    __m128d ca=_mm_set1_pd(cos(a));
    union dvec yz;
    yz.v=_mm_add_pd(_mm_mul_pd(ca,_mm_setr_pd(y(),z())),
			  _mm_mul_pd(sa,_mm_setr_pd(-z(),y())));
    xy.d[1]=yz.d[0];
    zx.d[0]=yz.d[1];
    return *this;
  } 

  // Rotate by angle a about the y-axis
  DVector3 RotateY(const double a){
    __m128d sa=_mm_set1_pd(sin(a));
    __m128d ca=_mm_set1_pd(cos(a));
    zx.v=_mm_add_pd(_mm_mul_pd(ca,zx.v),
		    _mm_mul_pd(sa,_mm_setr_pd(-x(),z())));
    xy.d[0]=zx.d[1];
    return *this;
  }

// Rotate by angle a about the axis specified by v.  The following ugly mess 
// of SIMD instructions does the following:
//            sa = sin(a), ca= cos(a)
//            dx=vx/|v|, dy=vy/|v|, dz=vz/|v|
// |x'| |ca+(1-ca)*dx*dx          (1-ca)*dx*dy-sa*dz    (1-ca)*dx*dz+sa*dy||x|
// |y'|=|   (1-ca)*dy*dx+sa*dz ca+(1-ca)*dy*dy          (1-ca)*dy*dz-sa*dx||y|
// |z'| |   (1-ca)*dz*dx-sa*dy    (1-ca)*dz*dy+sa*dx ca+(1-ca)*dz*dz      ||z|
//
  DVector3 Rotate(const double a,const DVector3 &v){
    if (a!=0){
      double vmag=v.Mag();
      if (vmag==0.0){
	//cerr << "axis vector has zero magnitude." <<endl;
      }
      else{
	__m128d xx=_mm_set1_pd(x());
	__m128d yy=_mm_set1_pd(y());
	__m128d zz=_mm_set1_pd(z());
	double sa=sin(a);
	double ca=cos(a);
	__m128d one_minus_ca=_mm_set1_pd(1.-ca);
	__m128d sa_zero=_mm_setr_pd(sa,0.);
	__m128d ca_zero=_mm_setr_pd(ca,0.);
	__m128d zero_sa=_mm_setr_pd(0.,sa);
	__m128d zero_ca=_mm_setr_pd(0.,ca);
	__m128d scale=_mm_set1_pd(1./vmag);
	__m128d dxdy=_mm_mul_pd(scale,v.GetVxy());
	__m128d dzdx=_mm_mul_pd(scale,v.GetVzx());
	__m128d dxdx=_mm_mul_pd(scale,_mm_set1_pd(v.x()));
	__m128d dydy=_mm_mul_pd(scale,_mm_set1_pd(v.y()));
	__m128d dzdz=_mm_mul_pd(scale,_mm_set1_pd(v.z()));
	__m128d AD=_mm_add_pd(_mm_add_pd(ca_zero,
					 _mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dxdy,dxdx))), 
			      _mm_mul_pd(zero_sa,dzdz));
	__m128d	BE=_mm_sub_pd(_mm_add_pd(zero_ca,
					 _mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dxdy,dydy))),
			      _mm_mul_pd(sa_zero,dzdz));
	__m128d CF=_mm_sub_pd(_mm_add_pd(_mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dxdy,dzdz)),
					 _mm_mul_pd(sa_zero,dydy)),
			      _mm_mul_pd(zero_sa,dxdx));
	__m128d IA=_mm_sub_pd(_mm_add_pd(zero_ca,
					 _mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dzdx,dxdx))),
			      _mm_mul_pd(sa_zero,dydy));
	__m128d	JB=_mm_sub_pd(_mm_add_pd(_mm_mul_pd(sa_zero,dxdx),
					 _mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dzdx,dydy))),
			      _mm_mul_pd(zero_sa,dzdz));
	__m128d KC=_mm_add_pd(_mm_add_pd(ca_zero,
					 _mm_mul_pd(one_minus_ca,
						    _mm_mul_pd(dzdx,dzdz))),
			      _mm_mul_pd(zero_sa,dydy));

	xy.v=_mm_add_pd(_mm_add_pd(_mm_mul_pd(AD,xx), _mm_mul_pd(BE,yy)),
			_mm_mul_pd(CF,zz));
	zx.v=_mm_add_pd(_mm_add_pd(_mm_mul_pd(IA,xx), _mm_mul_pd(JB,yy)),
			_mm_mul_pd(KC,zz));
      }
    }
    return *this;
  }
  
  // Take the dot product of "this" with the vector v.
  double Dot(const DVector3 &v) const{
    return (x()*v.x()+y()*v.y()+z()*v.z());
  };

  // Angle between two vectors
  double Angle(const DVector3 &v) const{
    double v1mag=Mag();
    double v2mag=v.Mag();
    if (v2mag>0. && v1mag > 0.){
      double costheta=Dot(v)/v1mag/v2mag;
      // Due to round-off errors, costheta could be epsilon greater than 1...
      if (costheta>1.) return 0.;
      if (costheta<-1.) return M_PI; 
      return acos(costheta);
    }
    return 0.;
  }

  // Print the components to the screen 
  void Print() const{
    cout << "DVector3 (x,y,z)=("<<setprecision(5)<<x() << ","
	 << setprecision(5) << y() << ","<<setprecision(5) << z() << "), " 
	 <<"(rho,theta,phi)=("<<setprecision(5) << Mag() << "," 
	 << setprecision(5) <<RAD2DEG*Theta() << "," 
	 << setprecision(5) << RAD2DEG*Phi() << ")" 
	  <<endl;   
  };


  // Routines to get the SIMD __mm128d vectors.
  __m128d GetVxy() const {return xy.v;};
  __m128d GetVzx() const {return zx.v;};

 private:  
  union dvec{
    __m128d v;
    double d[2];
  };
  union dvec xy;
  union dvec zx;
};

// Scale a vector by c
inline DVector3 operator*(const DVector3 &v1,const double c){
  __m128d scale=_mm_set1_pd(c);
  return DVector3(_mm_mul_pd(v1.GetVxy(),scale),_mm_mul_pd(v1.GetVzx(),scale));
};
// Scale a vector by c
inline DVector3 operator*(const double c,const DVector3 &v1){
  __m128d scale=_mm_set1_pd(c);
  return DVector3(_mm_mul_pd(v1.GetVxy(),scale),_mm_mul_pd(v1.GetVzx(),scale));
};


#endif // USE_SSE2
#endif // _DVector3_

