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
#include <align_16.h>
#define RAD2DEG 180./M_PI
using namespace std;

class DVector3{
 public:
  DVector3()
    : xy( ALIGNED_16_BLOCK_PTR(union dvec, 1, xy) )
    , zx(xy + 1)
    {
      xy->v=_mm_setzero_pd();
    zx->v=_mm_setzero_pd();
    };
  DVector3(double x, double y, double z)
  : xy( ALIGNED_16_BLOCK_PTR(union dvec, 1, xy) )
  , zx(xy + 1)
  {
    xy->v=_mm_setr_pd(x,y);
    zx->v=_mm_setr_pd(z,x); 
  };
  DVector3(__m128d v1, __m128d v2)
  : xy( ALIGNED_16_BLOCK_PTR(union dvec, 1, xy) )
  , zx(xy + 1)
  {
    xy->v=v1;
    zx->v=v2;
  }
    // Copy constructor
    DVector3(const DVector3 &v1): xy( ALIGNED_16_BLOCK_PTR(union dvec, 1, xy) )
    , zx(xy + 1){
      xy->v=v1.GetVxy();
      zx->v=v1.GetVzx();
    }

  ~DVector3(){};

  // Routines to set the components of the vector
  void SetXY(double x, double y){
    xy->v=_mm_setr_pd(x,y);
  }
  void SetXYZ(double x, double y, double z){
    xy->v=_mm_setr_pd(x,y);
    zx->v=_mm_setr_pd(z,x);
  }
  void SetMagThetaPhi(double p, double theta, double phi){
    double my_p=fabs(p);
    double pt=my_p*sin(theta); 
    xy->d[0]=zx->d[1]=pt*cos(phi);
    xy->d[1]=pt*sin(phi);
    zx->d[0]=my_p*cos(theta);
  }
  // Set phi keeping theta and magnitude fixed
  void SetPhi(double phi){
    double pt=Perp();
    xy->d[0]=zx->d[1]=pt*cos(phi);
    xy->d[1]=pt*sin(phi);
  }
  void SetX(double x){
    xy->d[0]=zx->d[1]=x;
  } 
  void SetY(double y){
    xy->d[1]=y;
  } 
  void SetZ(double z){
    zx->d[0]=z;
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
  
  double x() const {return xy->d[0];};
  double y() const {return xy->d[1];};
  double z() const {return zx->d[0];}; 
  double X() const {return xy->d[0];};
  double Y() const {return xy->d[1];};
  double Z() const {return zx->d[0];};
  double Px() const {return xy->d[0];};
  double Py() const {return xy->d[1];};
  double Pz() const {return zx->d[0];}; 

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
    return (xy->d[0]*xy->d[0]+xy->d[1]*xy->d[1]);
  }

  double Perp() const{
    return (sqrt(Perp2()));
  }
  double Pt() const{ return Perp();};

  double Mag2() const{
    return(xy->d[0]*xy->d[0]+xy->d[1]*xy->d[1]+zx->d[0]*zx->d[0]);
  }

  double Mag() const{
    return (sqrt(Mag2()));
  };

  // Cross product of "this"=(x,y,z) and v1
  //  |x'| |y*vz-z*vy|
  //  |y'|=|z*vx-x*vz|
  //  |z'| |x*vy-y*vx|
  DVector3 Cross(const DVector3 &v1) const {
    struct dvec1{
      __m128d v[2];
    };
    ALIGNED_16_BLOCK_WITH_PTR(struct dvec1, 1, p)
    struct dvec1 &yz=p[0];
    yz.v[0]=_mm_setr_pd(y(),z());
    yz.v[1]=_mm_setr_pd(v1.y(),v1.z());
    return DVector3(_mm_sub_pd(_mm_mul_pd(yz.v[0],v1.zx->v),
			       _mm_mul_pd(zx->v,yz.v[1])),
		    _mm_sub_pd(_mm_mul_pd(xy->v,yz.v[1]),
                               _mm_mul_pd(yz.v[0],v1.xy->v))  
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
  DVector3 SetMag(double c){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
    __m128d &scale=p[0];
    scale=_mm_set1_pd(c/Mag());
    xy->v=_mm_mul_pd(xy->v,scale);
    zx->v=_mm_mul_pd(zx->v,scale);
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
    xy->v=v1.xy->v; 
    zx->v=v1.zx->v;
    return *this;
  };

  // Addition
  DVector3 &operator+=(const DVector3 &v1){
    xy->v=_mm_add_pd(xy->v,v1.xy->v);
    zx->v=_mm_add_pd(zx->v,v1.zx->v);
    return *this;
  };
 
  // Subtraction
  DVector3 &operator-=(const DVector3 &v1){
    xy->v=_mm_sub_pd(xy->v,v1.xy->v);
    zx->v=_mm_sub_pd(zx->v,v1.zx->v);
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
  DVector3 &operator*=(double c){  
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
    __m128d &scale=p[0];
    scale=_mm_set1_pd(c);
    xy->v=_mm_mul_pd(GetVxy(),scale);
    zx->v=_mm_mul_pd(GetVzx(),scale);
    return *this;
  }
  
  // Unary minus
  DVector3 operator-() const{
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
    __m128d &zero=p[0];
    zero=_mm_set1_pd(0.);
    return DVector3(_mm_sub_pd(zero,xy->v),_mm_sub_pd(zero,zx->v));
  }

  // Rotate by angle a about the z-axis
  DVector3 RotateZ(double a){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 2, p)
    __m128d &ca=p[0];
    __m128d &sa=p[1];
    sa=_mm_set1_pd(sin(a));
    ca=_mm_set1_pd(cos(a));
    xy->v=_mm_add_pd(_mm_mul_pd(ca,xy->v),
		     _mm_mul_pd(sa,_mm_setr_pd(-y(),x())));
    zx->d[1]=x();
    return *this;
  } 

  // Rotate by angle a about the x-axis
  DVector3 RotateX(double a){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 2, p)
    __m128d &sa=p[0];
    __m128d &ca=p[1];
    sa=_mm_set1_pd(sin(a));
    ca=_mm_set1_pd(cos(a));
    union dvec2{
      __m128d v;
      double d[2];
    };
    ALIGNED_16_BLOCK_WITH_PTR(union dvec2, 1, yz)
    yz->v=_mm_add_pd(_mm_mul_pd(ca,_mm_setr_pd(y(),z())),
			  _mm_mul_pd(sa,_mm_setr_pd(-z(),y())));
    xy->d[1]=yz->d[0];
    zx->d[0]=yz->d[1];
    return *this;
  } 

  // Rotate by angle a about the y-axis
  DVector3 RotateY(double a){
    ALIGNED_16_BLOCK_WITH_PTR(__m128d, 2, p)
    __m128d &sa=p[0];
    __m128d &ca=p[1];
    sa=_mm_set1_pd(sin(a));
    ca=_mm_set1_pd(cos(a));
    zx->v=_mm_add_pd(_mm_mul_pd(ca,zx->v),
		    _mm_mul_pd(sa,_mm_setr_pd(-x(),z())));
    xy->d[0]=zx->d[1];
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
  DVector3 Rotate(double a,const DVector3 &v){
    if (a!=0){
      double vmag=v.Mag();
      if (vmag==0.0){
	//cerr << "axis vector has zero magnitude." <<endl;
      }
      else{
	double sa=sin(a);
	double ca=cos(a);
        ALIGNED_16_BLOCK_WITH_PTR(__m128d, 20, p)
        __m128d &xx=p[0];
        __m128d &yy=p[1];
        __m128d &zz=p[2];
        __m128d &one_minus_ca=p[3];
        __m128d &sa_zero=p[4];
        __m128d &ca_zero=p[5];
        __m128d &zero_sa=p[6];
        __m128d &zero_ca=p[7];
        __m128d &scale=p[8];
        __m128d &dxdy=p[9];
        __m128d &dzdx=p[10];
        __m128d &dxdx=p[11];
        __m128d &dydy=p[12];
        __m128d &dzdz=p[13];
        __m128d &AD=p[14];
        __m128d &BE=p[15];
        __m128d &CF=p[16];
        __m128d &IA=p[17];
        __m128d &JB=p[18];
        __m128d &KC=p[19];
	xx=_mm_set1_pd(x());
	yy=_mm_set1_pd(y());
	zz=_mm_set1_pd(z());
	one_minus_ca=_mm_set1_pd(1.-ca);
	sa_zero=_mm_setr_pd(sa,0.);
	ca_zero=_mm_setr_pd(ca,0.);
	zero_sa=_mm_setr_pd(0.,sa);
	zero_ca=_mm_setr_pd(0.,ca);
	scale=_mm_set1_pd(1./vmag);
	dxdy=_mm_mul_pd(scale,v.GetVxy());
	dzdx=_mm_mul_pd(scale,v.GetVzx());
	dxdx=_mm_mul_pd(scale,_mm_set1_pd(v.x()));
	dydy=_mm_mul_pd(scale,_mm_set1_pd(v.y()));
	dzdz=_mm_mul_pd(scale,_mm_set1_pd(v.z()));
	AD=_mm_add_pd(_mm_add_pd(ca_zero,
					_mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dxdy,dxdx))), 
			      _mm_mul_pd(zero_sa,dzdz));
	BE=_mm_sub_pd(_mm_add_pd(zero_ca,
					 _mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dxdy,dydy))),
			      _mm_mul_pd(sa_zero,dzdz));
	CF=_mm_sub_pd(_mm_add_pd(_mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dxdy,dzdz)),
					 _mm_mul_pd(sa_zero,dydy)),
			      _mm_mul_pd(zero_sa,dxdx));
	IA=_mm_sub_pd(_mm_add_pd(zero_ca,
					 _mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dzdx,dxdx))),
			      _mm_mul_pd(sa_zero,dydy));
	JB=_mm_sub_pd(_mm_add_pd(_mm_mul_pd(sa_zero,dxdx),
					 _mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dzdx,dydy))),
			      _mm_mul_pd(zero_sa,dzdz));
	KC=_mm_add_pd(_mm_add_pd(ca_zero,
					 _mm_mul_pd(one_minus_ca,
						_mm_mul_pd(dzdx,dzdz))),
			      _mm_mul_pd(zero_sa,dydy));

	xy->v=_mm_add_pd(_mm_add_pd(_mm_mul_pd(AD,xx),
			_mm_mul_pd(BE,yy)), _mm_mul_pd(CF,zz));
	zx->v=_mm_add_pd(_mm_add_pd(_mm_mul_pd(IA,xx),
			_mm_mul_pd(JB,yy)), _mm_mul_pd(KC,zz));
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
  __m128d GetVxy() const {return xy->v;};
  __m128d GetVzx() const {return zx->v;};

 private:  
  union dvec{
    __m128d v;
    double d[2];
  };
  ALIGNED_16_BLOCK(union dvec, 2, xy)
  union dvec* const zx;
};

// Scale a vector by c
inline DVector3 operator*(const DVector3 &v1, double c){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
  __m128d &scale=p[0];
  scale=_mm_set1_pd(c);
  return DVector3(_mm_mul_pd(v1.GetVxy(),scale),
                  _mm_mul_pd(v1.GetVzx(),scale));
}
// Scale a vector by c
inline DVector3 operator*(double c,const DVector3 &v1){
  ALIGNED_16_BLOCK_WITH_PTR(__m128d, 1, p)
  __m128d &scale=p[0];
  scale=_mm_set1_pd(c);
  return DVector3(_mm_mul_pd(v1.GetVxy(),scale),
                  _mm_mul_pd(v1.GetVzx(),scale));
}


#endif // USE_SSE2
#endif // _DVector3_

