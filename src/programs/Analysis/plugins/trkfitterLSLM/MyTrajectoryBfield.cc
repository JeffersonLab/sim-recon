#include <iostream>
#include <vector>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/Vector.h>

using namespace std;

#include "MyTrajectory.h"
#include "MyTrajectoryBfield.h"

#define TRACKING_RADIUS_MAX 60.0

MyTrajectoryBfield::MyTrajectoryBfield(const HepVector B_in,
				       int level)
  : MyTrajectory(level), nparams(5), params(nparams), BConst(B_in), delta(nparams, 0.0001), bfield(NULL),
    debug_level(level) {
  return;
}

MyTrajectoryBfield::MyTrajectoryBfield(const DMagneticFieldMap *bfield_in,
				       int level)
  : MyTrajectory(level), nparams(5), params(nparams), delta(nparams, 0.0001), bfield(bfield_in),
    debug_level(level) {
  return;
}

unsigned int MyTrajectoryBfield::getNumberOfParams() {
  return nparams;
}

void MyTrajectoryBfield::swim(const HepVector& param) {
  double xp0, z0, theta0, phi0, ptinv;
  xp0 = param(1);
  z0 = param(2);
  theta0 = param(3);
  phi0 = param(4);
  ptinv = param(5);
  MyTrajectoryBfield::swim(xp0, z0, theta0, phi0, ptinv);
  return;
}

void MyTrajectoryBfield::swim(double xp0, double z0, double theta,
			     double phi, double ptinv) {
  if (debug_level >= 3) cout << "starting B-field swim, xp0 = " << xp0
			     << " z0 = " << z0 << " theta = " << theta
			     << " phi = " << phi << " ptinv = " << ptinv << endl;
  // store input parameters as member data
  params(1) = xp0;
  params(2) = z0;
  params(3) = theta;
  params(4) = phi;
  params(5) = ptinv;
  checkClear();
  const double c = 29.9792458; // speed of light, cm/ns
  double m=0.139, dt=0.01; // gev/c^2, ns
  double cdt = c*dt;
  HepVector r_0(3), r_1(3), r_2(3), k(3), v_0(3);
  HepMatrix A(3,3,0), I(3,3,1), M(3,3), onePlusA(3,3);
  int ierr;
  double beta, k_1;
  if (ptinv != 0.0) {
    double p = 1.0/ptinv/sin(theta);
    double E = sqrt(m*m + p*p);
    double gamma = E/m;
    beta = sqrt(1.0 - 1.0/(gamma*gamma));
    double e;
    if (ptinv > 0.0) {e = 1.0;} else {e = -1.0;} // only allow singly
                                                  // charged particles
    k_1 = (e*1.60217653e-19)*(dt*1.0e-9)/(gamma*(m*1.78266181e-27));
    //    cout << "gamma = " << gamma << endl;
  } else {
    beta = 1.0;
    k_1 = 0.0;
  }
  //  cout << "beta = " << beta << " k_1 = " << k_1 << endl;
  double v = beta*c; // cm/ns
  double sinTheta = sin(theta);
  double cosTheta = cos(theta);
  double sinPhi = sin(phi);
  double cosPhi = cos(phi);
  v_0(1) = v*sinTheta*cosPhi;
  v_0(2) = v*sinTheta*sinPhi;
  v_0(3) = v*cosTheta;
  //  cout << "MyTrajectoryBfield::swim initial velocity:" << v_0;
  HepLorentzVector *thisVector;
  r_0(1) = xp0*sinPhi; // = xp0*cos(alpha) where alpha = phi - pi/2
  r_0(2) = -xp0*cosPhi; // = xp0*sin(alpha)
  r_0(3) = z0;
  thisVector = new HepLorentzVector(r_0(1), r_0(2), r_0(3), 0.0);
  traj.push_back(thisVector);
  r_1 = r_0 + v_0*dt;
  thisVector = new HepLorentzVector(r_1(1), r_1(2), r_1(3), cdt);
  traj.push_back(thisVector);
  double ctime = 2.0*cdt;
  HepVector B(3);
  for (int istep = 0; istep < 2000; istep++) {
    B = getField(r_1);
    if (debug_level > 2) cout << "MyTrajectoryBfield::swim B vector:" << B;
    k = 0.5*k_1*B;
    //cout << "MyTrajectoryBfield::swim k vector:" << k;
    A(1,2) = -k(3);
    A(1,3) = k(2);
    A(2,1) = k(3);
    A(2,3) = -k(1);
    A(3,1) = -k(2);
    A(3,2) = k(1);
    onePlusA = I + A;
    M = onePlusA.inverse(ierr)*(I - A); // can we take advantage of anti-sym?
    r_2 = r_1 + M*(r_1 - r_0);
    if (debug_level >=3 ) {
      if (istep%300 == 0) cout << istep << ' ' << r_2(1) << ' '
			       << r_2(2) << ' ' << r_2(3) << endl;
    }
    thisVector = new HepLorentzVector(r_2(1), r_2(2), r_2(3), ctime);
    traj.push_back(thisVector);
    r_0 = r_1; r_1 = r_2;
    ctime += cdt;
    if (sqrt(r_0(1)*r_0(1) + r_0(2)*r_0(2)) > TRACKING_RADIUS_MAX) break;
  }
  return;
}

HepVector MyTrajectoryBfield::getParams() {return params;}

HepVector MyTrajectoryBfield::getField(HepVector& r) {
  HepVector BThis(3);
  if (bfield) {
    bfield->GetField(r(1), r(2), r(3), BThis(1), BThis(2), BThis(3));
  } else {
    BThis = BConst;
  }
  return BThis;
}
