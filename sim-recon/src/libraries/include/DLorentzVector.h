// $Id$
//
//    File: DLorentzVector.h
//
// This header file is a replacement of a subset of TLorentzVector intended to 
// be consistent with the SIMDized version DVector3. 
//

#ifndef _DLorentzVector_
#define _DLorentzVector_

#ifndef USE_SSE2

#include <TLorentzVector.h>
typedef TLorentzVector DLorentzVector;

#else

#include "DVector3.h"
#include <math.h>
#include <emmintrin.h>
#include <iostream>
using namespace std;

class DLorentzVector{
 public:
  DLorentzVector(){
    mP.SetXYZ(0.,0.,0.);
    mE=0.;
  };
  DLorentzVector(const double x,const double y,const double z,const double t){
    mP.SetXYZ(x,y,z);
    mE=t;
  };
  DLorentzVector(const DVector3 &v, const double t){
    mP=v;
    mE=t;
  };
  ~DLorentzVector(){};
  void SetXYZT(const double x,const double y,const double z,const double t){
    mP.SetXYZ(x,y,z);
    mE=t;
  }
  // Set the 3-momentum or position part of the 4-vector
  void SetVect(const DVector3 &p){
    mP=p;
  }
  // Set the time or energy component
  void SetT(const double t){ mE=t;};
  // Set position components
  void SetX(const double x){mP.SetX(x);};
  void SetY(const double y){mP.SetY(y);};
  void SetZ(const double z){mP.SetZ(z);};
  
  // Routines to get position and time
  double X() const {return mP.x();};
  double Y() const {return mP.y();};
  double Z() const {return mP.z();}; 
  double T() const {return mE;};

  // Routine to get full 3-vector;
  DVector3 Vect() const {return mP;};

  // Routines to get momentum and energy
  double Px() const {return mP.x();};
  double Py() const {return mP.y();};
  double Pz() const {return mP.z();}; 
  double Pt() const {return mP.Perp();};
  double P() const {return mP.Mag();};
  double E() const {return mE;};
  double Energy() const {return mE;};

  // Spherical coordinates of spatial component
  double Rho() const { return mP.Mag();};

  // Angles
  double Theta() const {return mP.Theta();};
  double Phi() const {return mP.Phi();};

  // Kinematical quantities 
  double Beta() const { return P()/E();};
  double Mag2() const {return mE*mE-mP.Mag2();};
  double M() const{
    double mm = Mag2();
    return mm < 0.0 ? -sqrt(-mm) : sqrt(mm);
  }
  double M2() const {return Mag2();};
  double Mag() const {return M();};
  
  // Addition and subtraction
  DLorentzVector &operator+=(const DLorentzVector &v1){
    mP+=v1.Vect();
    mE+=v1.E();
    return *this;
  } 
  DLorentzVector &operator-=(const DLorentzVector &v1){
    mP-=v1.Vect();
    mE-=v1.E();
    return *this;
  }


  void Print() const{
    cout << "DLorentzVector (x,y,z,t)=(" << X() << "," << Y() << "," << Z()
	 << "," << T() << ")" << endl;

  };

 private:
  DVector3 mP;  // momentum or position vector
  double mE;  // Energy or time component
};

// Addition 
inline DLorentzVector operator+(const DLorentzVector &v1,const DLorentzVector &v2){
  return DLorentzVector(v1.Vect()+v2.Vect(),v1.E()+v2.E());
}
//Subtraction 
inline DLorentzVector operator-(const DLorentzVector &v1,const DLorentzVector &v2){
  return DLorentzVector(v1.Vect()-v2.Vect(),v1.E()-v2.E());
}

#endif // USE_SSE2

#endif // _DLorentzVector_

