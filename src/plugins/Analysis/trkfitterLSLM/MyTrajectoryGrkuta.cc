#include "MyTrajectoryGrkuta.h"
#define M_CHARGED_PI 0.13957018
#define c 29.9792548

#define TRACKING_RADIUS_MAX 60.0
#define TIME_STEP_SIZE 0.03

MyTrajectoryGrkuta* globalGrkutaPtr;

MyTrajectoryGrkuta::MyTrajectoryGrkuta(const HepVector B_in, int level) : MyTrajectoryBfield(B_in), debug_level(level) {
}

MyTrajectoryGrkuta::MyTrajectoryGrkuta(const DMagneticFieldMap *bfield_in,
				       int level)
  : MyTrajectoryBfield(bfield_in, level), debug_level(level) {
  return;
}

extern "C" {
  void grkuta_(float* charge, float* step, float* vect, float* vout);
  void gufld_(float* x, float* f);
}

// The following member function is a cut and paste from
// MyTrajectoryBfield. Is there a better way?

void MyTrajectoryGrkuta::swim(const HepVector& param) {
  double xp0, z0, theta0, phi0, ptinv;
  xp0 = param(1);
  z0 = param(2);
  theta0 = param(3);
  phi0 = param(4);
  ptinv = param(5);
  MyTrajectoryGrkuta::swim(xp0, z0, theta0, phi0, ptinv);
  return;
}

void MyTrajectoryGrkuta::swim(double xp0, double z0, double theta, double phi,
			 double ptinv) {

  // store input parameters as member data
  params(1) = xp0;
  params(2) = z0;
  params(3) = theta;
  params(4) = phi;
  params(5) = ptinv;
  checkClear();

  double sinTheta = sin(theta);
  double cosTheta = cos(theta);
  double sinPhi = sin(phi);
  double cosPhi = cos(phi);
  double xStart = xp0*sinPhi; // = xp0*cos(alpha) where alpha = phi - pi/2
  double yStart = -xp0*cosPhi; // = xp0*sin(alpha)
  double zStart = z0;

  HepLorentzVector* thisVectorPtr;
  float charge, step, vect[7], vout[7];
  if (ptinv > 0.0) {charge = 1.0;} else {charge = -1.0;}
  vect[0] = xStart; vect[1] = yStart; vect[2] = zStart;
  vect[3] = sinTheta*cosPhi; vect[4] = sinTheta*sinPhi; vect[5] = cosTheta;
  double ptot = abs(1.0/(ptinv*sinTheta));
  vect[6] = ptot;
  double energy = sqrt(M_CHARGED_PI*M_CHARGED_PI + ptot*ptot);
  double gamma = energy/M_CHARGED_PI;
  double beta = sqrt(1.0 - 1.0/(gamma*gamma));
  double tStep = TIME_STEP_SIZE;
  double ctStep = c*tStep;
  step = (float)(beta*ctStep);
  double ct = 0.0;
  globalGrkutaPtr = this;
  thisVectorPtr = new HepLorentzVector;
  thisVectorPtr->setX(vect[0]);
  thisVectorPtr->setY(vect[1]);
  thisVectorPtr->setZ(vect[2]);
  thisVectorPtr->setT(ct);
  traj.push_back(thisVectorPtr);
  for (int i = 0; i < 2000; i++) {
    grkuta_(&charge, &step, vect, vout);
    ct += ctStep;
    thisVectorPtr = new HepLorentzVector;
    thisVectorPtr->setX(vout[0]);
    thisVectorPtr->setY(vout[1]);
    thisVectorPtr->setZ(vout[2]);
    thisVectorPtr->setT(ct);
    if (debug_level > 3) cout << setprecision(14) << "MyTrajectoryGrcuta::swim: i = " << i << " LorentzVector = " << *thisVectorPtr << endl;
    traj.push_back(thisVectorPtr);
    if (sqrt(vout[0]*vout[0] + vout[1]*vout[1]) > TRACKING_RADIUS_MAX) break;
    for (int j = 0; j < 7; j++) {vect[j] = vout[j];}
  }
  return;
}

void gufld_(float* x, float* f) {
  const double kGperT = 10.0; // kilogauss per Tesla
  HepVector r(3);
  r(1) = (double)x[0];
  r(2) = (double)x[1];
  r(3) = (double)x[2];
  HepVector B = globalGrkutaPtr->getField(r);
  f[0] = kGperT*B(1);
  f[1] = kGperT*B(2);
  f[2] = kGperT*B(3);
  return;
}
